namespace imajuscule::fft {
    
    struct CallCost {
        using MeasureCall = std::function<profiling::CpuDuration(int64_t, int64_t, double &)>;
        
        CallCost(MeasureCall m)
        : measure_call(m)
        {}
        
        double operator ()(int64_t const sz) {
            if(!sz) {
                return 0.;
            }
            std::lock_guard<std::mutex> l(mut);
            
            auto it = stats.find(sz);
            if(it != stats.end()) {
                return it->second;
            }
            double const t = forceTiming(sz);
            stats[sz] = t;
            return t;
        }

        double forceTiming(int64_t const sz) const
        {
            constexpr int minMicroSeconds = 4;
            int64_t ntests = 1;
            // We have a timing granularity of 1 micro second, so we double the number of tests
            // until the timing is well-bigger than the granularity
            std::optional<double> minT;
            double sideEffect = 0.;
start_again:
            for(int i=0; i<4; ++i) {
                auto microseconds = measure_call(sz, ntests, sideEffect).count();
                if(microseconds < minMicroSeconds) {
                    ntests *= 2;
                    minT.reset();
                    goto start_again;
                }
                double localTPerTest = microseconds / (1000000. * static_cast<double>(ntests));
                if(i==0) {
                    // don't use the timing, it was just to warm up the cache
                    continue;
                }
                minT = minT ? std::min(*minT,localTPerTest) : localTPerTest;
            }
            return *minT + std::min(std::abs(sideEffect), 0.000000000001);
        }
        
    private:
        std::unordered_map<int, double> stats;
        std::mutex mut;
        
        MeasureCall measure_call;
    };
    
    template<typename Tag, typename T>
    struct RealSignalCosts {
        using Impl = RealSignal_<Tag, T>;
        using type = typename Impl::type;
        using value_type = typename type::value_type;

        inline static CallCost cost_add_scalar_multiply { [](int64_t sz, int64_t ntests, double & sideEffect) {
            std::vector<std::array<type, 3>> va;
            va.resize(ntests);
            double d = 1.;
            for(auto & a:va)
            {
                for(auto & v:a) {
                    v.resize(sz,
                             value_type{static_cast<T>(d)});
                    ++d;
                }
            }
            
            using namespace profiling;
            CachePolluter flushCpuCaches(sideEffect);
            auto duration = measure_thread_cpu_one([&va, sz](){
                for(auto & a:va) {
                    Impl::add_scalar_multiply(a[0].begin(),
                                              a[1].begin(),
                                              a[2].begin(),
                                              0.7,
                                              sz);
                }
            });
            for(auto & a:va)
            {
                using std::abs; // so that abs can be std::abs or imajuscule::abs (for complex)
                sideEffect += abs(std::accumulate(a[0].begin(),
                                                  a[0].end(),
                                                  value_type{})) / static_cast<T>(a[0].size());
            }
            return duration;
        }};
        
        inline static CallCost cost_copy { [](int64_t sz, int64_t ntests, double & sideEffect) {
            std::vector<std::array<type, 2>> va;
            va.resize(ntests);
            {
                int i=0;
                for(auto & a : va)
                {
                    for(auto & v:a) {
                        ++i;
                        v.resize(sz, value_type{static_cast<T>(i)});
                    }
                }
            }
            using namespace profiling;
            CachePolluter flushCpuCaches(sideEffect);
            auto duration = measure_thread_cpu_one([&va, sz](){
                for(auto & a : va) {
                    Impl::copy(a[0].begin(),
                               a[1].begin(),
                               sz);
                }
            });
            using std::abs;
            for(auto & a : va)
            {
                sideEffect += abs(std::accumulate(a[0].begin(),
                                                  a[0].end(),
                                                  value_type{})) / static_cast<T>(a[0].size());
            }
            return duration;
        }};

        inline static CallCost cost_zero_n_raw { [](int64_t sz, int64_t ntests, double & sideEffect) {
            std::vector<type> va;
            va.resize(ntests);
            for(auto & a:va) {
                a.resize(sz, value_type(1.));
            }
            
            for(auto & a:va) {
                using std::abs;
                sideEffect += abs(a[0]);
            }
            
            using namespace profiling;
            CachePolluter flushCpuCaches(sideEffect);
            auto duration = measure_thread_cpu_one([&va, sz](){
                for(auto & a : va) {
                    Impl::zero_n_raw(&a[0], sz);
                }
            });
            
            for(auto & a:va) {
                using std::abs;
                sideEffect += abs(a[sz-1]);
            }
            return duration;
        }};
        
    };
    

    template<typename Tag, typename T>
    struct RealFBinsCosts {
        // we use a simple allocator to compute costs, one that doesn't take too much space
        // (hence, not page-aligned allocator) and that is simple to use (hence, not monotonic allocator)
        using Impl = RealFBins_<Tag, T, a64::Alloc>;

        using type = typename Impl::type;
        using value_type = typename type::value_type;

        static constexpr auto getFirstReal = Impl::getFirstReal;
        static constexpr auto setFirstReal = Impl::setFirstReal;

        inline static CallCost cost_mult_assign { [](int64_t sz, int64_t ntests, double & sideEffect) {
            std::vector<std::array<type, 2>> va;
            va.resize(ntests);
            int i=0;
            for(auto & a:va) {
                for(auto & v:a) {
                    ++i;
                    v.resize(sz);
                    setFirstReal(v, static_cast<T>(i));
                }
            }
            using namespace profiling;
            CachePolluter flushCpuCaches(sideEffect);
            auto duration = measure_thread_cpu_one([&va, sz](){
                for(auto & a:va) {
                    Impl::mult_assign(a[0],
                                      a[1]);
                }
            });
            for(auto & a:va) {
                sideEffect += getFirstReal(a[0]);
            }
            return duration;
        }};
        
        inline static CallCost cost_multiply { [](int64_t sz, int64_t ntests, double & sideEffect) {
            std::vector<std::array<type, 3>> va;
            va.resize(ntests);
            for(auto & a:va) {
                int i=0;
                for(auto & v:a) {
                    ++i;
                    v.resize(sz);
                    setFirstReal(v, static_cast<T>(i));
                }
            }
            using namespace profiling;
            auto duration = measure_thread_cpu_one([&va, sz](){
                for(auto & a:va) {
                    Impl::multiply(a[0],
                                   a[1],
                                   a[2]);
                }
            });
            for(auto & a:va) {
                sideEffect += getFirstReal(a[0]);
            }
            return duration;
        }};
        
        inline static CallCost cost_multiply_add { [](int64_t sz, int64_t ntests, double & sideEffect) {
            std::vector<std::array<type, 3>> va;
            va.resize(ntests);
            {
                int i=0;
                for(auto & a:va) {
                    for(auto & v:a) {
                        ++i;
                        v.resize(sz);
                        setFirstReal(v, static_cast<T>(i));
                    }
                }
            }
            
            using namespace profiling;
            CachePolluter flushCpuCaches(sideEffect);
            auto duration = measure_thread_cpu_one([&va](){
                for(auto & a:va) {
                    Impl::multiply_add(a[0],
                                       a[1],
                                       a[2]);
                }
            });

            for(auto & a:va) {
                sideEffect += getFirstReal(a[0]);
            }
            return duration;
        }};
    };
    
    template<typename Tag, typename T>
    struct AlgoCosts {
        using RealInput = typename RealSignal_ <Tag, T>::type;
        
        // we use a simple allocator to compute costs, one that doesn't take too much space
        // (hence, not page-aligned allocator) and that is simple to use (hence, not monotonic allocator)
        template<typename TT>
        using Allocator = a64::Alloc<TT>;
        
        using RealFBins = RealFBins_<Tag, T, Allocator>;
        using RealFBins_t  = typename RealFBins::type;
        using Contexts = fft::Contexts_<Tag, T>;
        using Algo = Algo_<Tag, T>;

        static constexpr auto setFirstReal = RealFBins::setFirstReal;
        static constexpr auto getFirstReal = RealFBins::getFirstReal;

        inline static CallCost cost_fft_inverse { [](int64_t sz, int64_t ntests, double & sideEffect) {
            Assert(sz);
            Algo a(Contexts::getInstance().getBySize(sz));

            std::vector<RealInput> vinput;
            vinput.resize(ntests);
            for(auto &input:vinput) {
                input.resize(sz);
            }

            std::vector<RealFBins_t> vf;
            vf.resize(ntests);
            {
                int i=1;
                for(auto & f:vf) {
                    f.resize(sz);
                    setFirstReal(f, static_cast<T>(i));
                    ++i;
                }
            }

            using namespace profiling;
            CachePolluter flushCpuCaches(sideEffect);
            auto duration = measure_thread_cpu_one([&a, &vf, &vinput, sz, ntests](){
                for(int i=0; i<ntests; ++i) {
                    a.inverse(vf[i].data(), vinput[i], sz);
                }
            });
            
            for(auto &input:vinput) {
                using std::abs;
                sideEffect += abs(input[0]);
            }
            return duration;
        }};
                    
        inline static CallCost cost_fft_forward { [](int64_t sz, int64_t ntests, double & sideEffect) {
            Assert(sz);
            Algo a(Contexts::getInstance().getBySize(sz));

            std::vector<RealFBins_t> vf;
            vf.resize(ntests);
            for(auto & f:vf) {
                f.resize(sz);
            }
            
            std::vector<RealInput> vinput;
            vinput.resize(ntests);
            {
                int i=1;
                for(auto & input:vinput) {
                    input.resize(sz);
                    input[0] = typename RealInput::value_type(static_cast<T>(i));
                    ++i;
                }
            }

            using namespace profiling;
            CachePolluter flushCpuCaches(sideEffect);
            auto duration = measure_thread_cpu_one([&a, &vf, &vinput, sz, ntests](){
                for(int i=0; i<ntests; ++i) {
                    a.forward(vinput[i].begin(), vf[i].data(), sz);
                }
            });
            
            for(auto &f:vf) {
                sideEffect += getFirstReal(f);
            }
            return duration;
        }};
    };
} // NS imajuscule::fft

namespace imajuscule {

    struct XFFtCostFactors {
        XFFtCostFactors(std::map<int, float> const & multiplicators = {})
        : multiplicators(multiplicators)
        {}
        
        float getCostMultiplicator(int fftSize) const {
            Assert(!fftSize || is_power_of_two(fftSize));
            auto it = multiplicators.find(fftSize/2);
            if(it == multiplicators.end()) {
                return 1.f;
            }
            return it->second;
        }
        
        void setMultiplicator(int sz, float factor) {
            multiplicators[sz] = factor;
        }
    private:
        // fft half-size -> multiplicator
        std::map<int, float> multiplicators;
    };
} // NS imajuscule

