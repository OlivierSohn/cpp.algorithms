namespace imajuscule::fft {
    
    template<typename Tag, typename T>
    struct RealSignalCosts {
        using Impl = RealSignal_<Tag, T>;
        using type = typename Impl::type;
        using value_type = typename type::value_type;

        static double cost_add_scalar_multiply(int64_t sz) {
            static std::unordered_map<int, double> stats;
            static std::mutex mut;
            std::lock_guard<std::mutex> l(mut);
            
            auto it = stats.find(sz);
            if(it != stats.end()) {
                return it->second;
            }
            double t = nocache_cost_add_scalar_multiply(sz);
            stats[sz] = t;
            return t;
        }
        static double nocache_cost_add_scalar_multiply(int64_t sz) {
            std::array<type, 3> a;
            for(auto & v:a) {
                v.resize(sz, value_type{1.});
            }
            int64_t ntests = std::max(static_cast<int64_t>(1),
                                      10000 / sz);
            using namespace profiling;
            auto duration = measure_thread_cpu_one([&a, sz, ntests](){
                std::array<int, 3> idx {0,1,2};
                
                for(int i=0; i<ntests; i++) {
                    Impl::add_scalar_multiply(a[idx[0]].begin(),
                                              a[idx[1]].begin(),
                                              a[idx[2]].begin(),
                                              0.7,
                                              sz);
                    std::rotate(idx.begin(),
                                idx.begin()+1,
                                idx.end());
                }
            });
            using std::abs; // so that abs can be std::abs or imajuscule::abs (for complex)
            T avg = abs(std::accumulate(a[2].begin(), a[2].end(), value_type{})) / static_cast<T>(a[2].size());
            return duration.count() / (1000000. * static_cast<double>(ntests) + std::min(static_cast<T>(0.000001), avg));
        }
        

        static double cost_copy(int64_t sz) {
            static std::unordered_map<int, double> stats;
            static std::mutex mut;
            std::lock_guard<std::mutex> l(mut);
            
            auto it = stats.find(sz);
            if(it != stats.end()) {
                return it->second;
            }
            double t = nocache_cost_copy(sz);
            stats[sz] = t;
            return t;
        }
        static auto nocache_cost_copy(int64_t sz) {
            std::array<type, 2> a;
            {
                int i=0;
                for(auto & v:a) {
                    ++i;
                    v.resize(sz, value_type{static_cast<T>(i)});
                }
            }
            int64_t ntests = std::max(static_cast<int64_t>(1),
                                      10000 / sz);
            using namespace profiling;
            auto duration = measure_thread_cpu_one([&a, sz, ntests](){
                std::array<int, 3> idx {0,1};
                for(int i=0; i<ntests; i++) {
                    Impl::copy(a[idx[0]].begin(),
                               a[idx[1]].begin(),
                               sz);
                    std::rotate(idx.begin(),
                                idx.begin()+1,
                                idx.end());
                }
            });
            using std::abs;
            T avg = abs(std::accumulate(a[1].begin(), a[1].end(), value_type{})) / static_cast<T>(a[1].size());
            return duration.count() / (1000000. * static_cast<double>(ntests) + std::min(static_cast<T>(0.000001), avg));
        }

        static double cost_zero_n_raw(int64_t sz) {
            static std::unordered_map<int, double> stats;
            static std::mutex mut;
            std::lock_guard<std::mutex> l(mut);
            
            auto it = stats.find(sz);
            if(it != stats.end()) {
                return it->second;
            }
            double t = nocache_cost_zero_n_raw(sz);
            stats[sz] = t;
            return t;
        }
        static auto nocache_cost_zero_n_raw(int64_t sz) {
            type a;
            a.resize(sz, value_type(1.));
            
            int64_t ntests = std::max(static_cast<int64_t>(1),
                                      10000 / sz);
            using namespace profiling;
            T sum = 0;
            auto duration = measure_thread_cpu_one([&sum, &a, sz, ntests](){
                for(int i=0; i<ntests; i++) {
                    Impl::zero_n_raw(&a[0], sz);
                    if(sz) {
                        a[sz-1] += value_type(1.);
                        using std::abs;
                        sum += abs(a[sz-1]);
                    }
                }
            });
            sum = std::abs(sum);
            return duration.count() / (1000000. * static_cast<double>(ntests) + std::min(static_cast<T>(0.000001), sum));
        }
        
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

        static double cost_mult_assign(int64_t sz) {
            static std::unordered_map<int, double> stats;
            static std::mutex mut;
            std::lock_guard<std::mutex> l(mut);
            
            auto it = stats.find(sz);
            if(it != stats.end()) {
                return it->second;
            }
            double t = nocache_cost_mult_assign(sz);
            stats[sz] = t;
            return t;
        }
        static auto nocache_cost_mult_assign(int64_t sz) {
            std::array<type, 2> a;
            int i=0;
            for(auto & v:a) {
                ++i;
                v.resize(sz);
                setFirstReal(v, static_cast<T>(i));
            }
            int64_t ntests = std::max(static_cast<int64_t>(1),
                                      10000 / sz);
            using namespace profiling;
            auto duration = measure_thread_cpu_one([&a, sz, ntests](){
                for(int i=0; i<ntests; ++i) {
                    Impl::mult_assign(a[0],
                                      a[1]);
                }
            });
            T sum = std::abs(getFirstReal(a[0]));
            return duration.count() / (1000000. * static_cast<double>(ntests) + std::min(static_cast<T>(0.000001), sum));
        }
        
        static double cost_multiply(int64_t sz) {
            static std::unordered_map<int, double> stats;
            static std::mutex mut;
            std::lock_guard<std::mutex> l(mut);
            
            auto it = stats.find(sz);
            if(it != stats.end()) {
                return it->second;
            }
            double t = nocache_cost_multiply(sz);
            stats[sz] = t;
            return t;
        }
        static auto nocache_cost_multiply(int64_t sz) {
            std::array<type, 3> a;
            int i=0;
            for(auto & v:a) {
                ++i;
                v.resize(sz);
                setFirstReal(v, static_cast<T>(i));
            }
            int64_t ntests = std::max(static_cast<int64_t>(1),
                                      10000 / sz);
            using namespace profiling;
            auto duration = measure_thread_cpu_one([&a, sz, ntests](){
                std::array<int, 3> idx {0,1,2};
                for(int i=0; i<ntests; ++i) {
                    Impl::multiply(a[idx[0]],
                                   a[idx[1]],
                                   a[idx[2]]);

                    std::rotate(idx.begin(),
                                idx.begin()+1,
                                idx.end());
                }
            });
            T sum {};
            for(auto const & aa:a) {
                sum += getFirstReal(aa);
            }
            return duration.count() / (1000000. * static_cast<double>(ntests) + std::min(static_cast<T>(0.000001), std::abs(sum)));
        }
        
        static double cost_multiply_add(int64_t sz) {
            static std::unordered_map<int, double> stats;
            static std::mutex mut;
            std::lock_guard<std::mutex> l(mut);
            
            auto it = stats.find(sz);
            if(it != stats.end()) {
                return it->second;
            }
            double t = nocache_cost_multiply_add(sz);
            stats[sz] = t;
            return t;
        }
        static auto nocache_cost_multiply_add(int64_t sz) {
            std::array<type, 3> a;
            int i=0;
            for(auto & v:a) {
                ++i;
                v.resize(sz);
                setFirstReal(v, static_cast<T>(i));
            }
            int64_t ntests = std::max(static_cast<int64_t>(1),
                                      10000 / sz);
            using namespace profiling;
            auto duration = measure_thread_cpu_one([&a, sz, ntests](){
                std::array<int, 3> idx {0,1,2};
                for(int i=0; i<ntests; ++i) {
                    Impl::multiply_add(a[idx[0]],
                                       a[idx[1]],
                                       a[idx[2]]);

                    std::rotate(idx.begin(),
                                idx.begin()+1,
                                idx.end());
                }
            });
            T sum {};
            for(auto const & aa:a) {
                sum += getFirstReal(aa);
            }
            return duration.count() / (1000000. * static_cast<double>(ntests) + std::min(static_cast<T>(0.000001), std::abs(sum)));
        }
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

        static double cost_fft_inverse(int64_t sz) {
            static std::unordered_map<int, double> stats;
            static std::mutex mut;
            std::lock_guard<std::mutex> l(mut);
            
            auto it = stats.find(sz);
            if(it != stats.end()) {
                return it->second;
            }
            double t = nocache_cost_fft_inverse(sz);
            stats[sz] = t;
            return t;
        }
        static auto nocache_cost_fft_inverse(int64_t sz) {
            if(sz == 0) {
                return 0.;
            }
            RealFBins_t f;
            RealInput input;
            f.resize(sz);
            input.resize(sz);
            Algo a(Contexts::getInstance().getBySize(sz));

            int64_t ntests = std::max(static_cast<int64_t>(1),
                                      10000 / sz);

            using namespace profiling;
            T sum = 0;
            auto duration = measure_thread_cpu_one([&sum, &a, &f, &input, sz, ntests](){
                for(int i=0; i<ntests; ++i) {
                    setFirstReal(f, static_cast<T>(i));
                    a.inverse(f.data(), input, sz);
                    using std::abs;
                    sum += abs(input[0]);
                }
            });
            sum = std::abs(sum);
            return duration.count() / (1000000. * static_cast<double>(ntests) + std::min(static_cast<T>(0.000001), sum));
        }
                    
        static double cost_fft_forward(int64_t sz) {
            static std::unordered_map<int, double> stats;
            static std::mutex mut;
            std::lock_guard<std::mutex> l(mut);
            
            auto it = stats.find(sz);
            if(it != stats.end()) {
                return it->second;
            }
            double t = nocache_cost_fft_forward(sz);
            stats[sz] = t;
            return t;
        }
        static auto nocache_cost_fft_forward(int64_t sz) {
            if(sz == 0) {
                return 0.;
            }
            RealFBins_t f;
            RealInput input;
            f.resize(sz);
            input.resize(sz);
            Algo a(Contexts::getInstance().getBySize(sz));
            
            int64_t ntests = std::max(static_cast<int64_t>(1),
                                      10000 / sz);

            using namespace profiling;
            T sum = 0;
            auto duration = measure_thread_cpu_one([&sum, &a, &f, &input, sz, ntests](){
                for(int i=0; i<ntests; ++i) {
                    input[0] = typename decltype(input)::value_type(static_cast<T>(i));
                    a.forward(input.begin(), f.data(), sz);
                    sum += getFirstReal(f);
                }
            });
            sum = std::abs(sum);
            return duration.count() / (1000000. * static_cast<double>(ntests) + std::min(static_cast<T>(0.000001), sum));
        }
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

