namespace imajuscule::fft {

    template<typename M, typename F>
    void forwardPass(M & stats, F f) {
        auto it1 = stats.begin();
        if(it1 == stats.end()) {
            return;
        }
        
        for(auto it2 = std::next(it1);
            it2 != stats.end();
            ++it1, ++it2) {
            f(it1, it2);
        }
    }
    template<typename M, typename F>
    void backwardPass(M & stats, F f) {
        auto it1 = stats.rbegin();
        if(it1 == stats.rend()) {
            return;
        }
        
        for(auto it2 = std::next(it1);
            it2 != stats.rend();
            ++it1, ++it2) {
            f(it2, it1);
        }
    }

    struct CallCost {
        using MeasureCall = std::function<profiling::CpuDuration(int64_t, int64_t, double &)>;
        
        static constexpr int maxAggregateSize = 10;

        CallCost(MeasureCall m)
        : measure_call(m)
        {}
        
        double operator ()(int64_t const sz) {
            if(!sz) {
                return 0.;
            }
            std::lock_guard<std::mutex> l(mut);
            
            addMeasures(sz);
            return stats[sz].average();
        }
        
        struct Measure {
            double time_by_test;
            int ntests;
        };
        
        struct AggregatedMeasure {
            void include(Measure m) {
                measures.push_back(m);
                avg = computeAvg();
            }
            double average() const {
                return avg;
            }
            int size() const {
                return measures.size();
            }
        private:
            double avg;
            std::vector<Measure> measures;
            double computeAvg() const {
                double avg = 0.;
                for(auto const & m:measures)
                {
                    avg += m.time_by_test;
                }
                return avg/measures.size();
            }
        };

        template<typename F>
        void forEachCost(F f) const {
            std::lock_guard<std::mutex> l(mut);
            for(auto const & [sz, cost] : stats) {
                f(sz, cost);
            }
        }
        
        void clear() {
            std::lock_guard<std::mutex> l(mut);
            stats.clear();
        }
        
    private:
        std::map<int, AggregatedMeasure> stats;
        mutable std::mutex mut;
        
        MeasureCall measure_call;

        /*
         if the user asked for sz 8, we will automatically add sz 4, 8, 16.
         
         and if sz 128 and 256 are already existing,
         we will also add sz 32 and 64.
         
         Then we will (for some iterations, see maxAggregateSize)
         redo the measurements so as to try to enforce the monotonicity constraint.
         */
        void addMeasures(int const sz) {
            Assert(sz);
            if(0==stats.count(sz)) {
                stats[sz].include(forceTiming(sz));
            }
            auto it = stats.upper_bound(sz);
            
            // We want to have at least one measure every power of 2 on contiguous powers of 2.
            
            // We fill the gaps upward
            {
                int nextSz = ceil_power_of_two(sz+1);
                if(it != stats.end()) {
                    int const nextSzExisting = it->first;
                    while(nextSz < nextSzExisting) {
                        stats[nextSz].include(forceTiming(nextSz));
                        nextSz *= 2;
                    }
                }
                else {
                    stats[nextSz].include(forceTiming(nextSz));
                }
            }
            // We fill the gaps downward
            {
                int prevSz = floor_power_of_two(sz-1);
                if(!stats.empty() && it != stats.begin()) {
                    --it;
                    int const prevSzExisting = it->first;
                    while(prevSz && (prevSz > prevSzExisting)) {
                        stats[prevSz].include(forceTiming(prevSz));
                        prevSz /= 2;
                    }
                }
                else if(prevSz) {
                    stats[prevSz].include(forceTiming(prevSz));
                }
            }
            
            enforceMonotonicConstraint();
        }
        
        void enforceMonotonicConstraint()
        {
            while(true)
            {
                bool stable = true;
                backwardPass(stats,[this, &stable](auto it1, auto it2){
                    if(it1->second.size() >= maxAggregateSize) {
                        return;
                    }
                    if(it1->second.average() > it2->second.average()) {
                        // assume the measure of higher size is right
                        it1->second.include(forceTiming(it1->first));
                        stable = false;
                    }
                });
                forwardPass(stats,[this, &stable](auto it1, auto it2){
                    if(it2->second.size() >= maxAggregateSize) {
                        return;
                    }
                    if(it1->second.average() > it2->second.average()) {
                        // assume the measure of lower size is right
                        it2->second.include(forceTiming(it2->first));
                        stable = false;
                    }
                });
                if(stable) {
                    return;
                }
            }
        }

        Measure forceTiming(int64_t const sz) const
        {
            constexpr int minMicroSeconds = 2;
            int ntests = 1;
            // We have a timing granularity of 1 micro second, so we double the number of tests
            // until the timing is well-bigger than the granularity
            std::optional<double> minT;
            double sideEffect = 0.;
            constexpr int nDropped = 1; // to warm-up the cache
            constexpr int nUsed = 3;
            constexpr int nTests = nDropped + nUsed;
            
        start_again:

            for(int i=0; i<nTests; ++i) {
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
            return {
                *minT + std::min(std::abs(sideEffect), 0.000000000001),
                ntests
            };
        }
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

