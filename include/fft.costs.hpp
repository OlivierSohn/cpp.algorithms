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
        using MeasureCall = std::function<profiling::CpuDuration(std::function<void(double&)>, int64_t, int64_t, double &)>;

    private:
        /*
         The count of measurements is in the interval :
         [nConstecutiveMeasures, nConstecutiveMeasures * maxAggregateSize]
         */
#ifdef NDEBUG
        static constexpr int nConstecutiveMeasures = 3;
        static constexpr int maxAggregateSize = 4;
#else
        static constexpr int nConstecutiveMeasures = 1;
        static constexpr int maxAggregateSize = 1;
#endif
        

    public:
        
        CallCost(MeasureCall m)
        : measure_call(m)
        {}
        
        double operator ()(int64_t const sz) {
            if(!sz) {
                return 0.;
            }
            std::lock_guard<std::mutex> l(mut);
            
            if(is_power_of_two(sz)) {
                int const nAdded = addMeasures(sz);
                if(nAdded) {
                    enforceMonotonicConstraint();
                }
                return stats[sz].minimum();
            }
            else {
                // do a linear interpolation
                
                int const ceil2 = ceil_power_of_two(sz);
                int const floor2 = floor_power_of_two(sz);
                int const middle2 = (ceil2 + floor2) / 2;
                Assert(ceil2 > middle2 && middle2 > floor2);
                Assert(ceil2 > sz && sz > floor2);
                int const nAdded = addMeasures(middle2) + addMeasures(floor2) + addMeasures(ceil2);
                if(nAdded) {
                    enforceMonotonicConstraint();
                }
                double const ceilVal = stats[ceil2].minimum();
                double const floorVal = stats[floor2].minimum();
                double const midValMinimumBound = (ceilVal + floorVal) / 2.; // so that interpolation > floorVal
                double const midVal = std::max(midValMinimumBound,
                                               stats[middle2].minimum());
                double const interpolation = midVal + (ceilVal - midVal) * (sz-middle2)/static_cast<double>(ceil2-middle2);
                Assert(ceilVal >= interpolation && interpolation >= floorVal);
                return interpolation;
            }
        }
        
        struct Measure {
            double time_by_test;
            int ntests;
        };
        
        struct AggregatedMeasure {
            void include(Measure m) {
                Assert(size() < maxAggregateSize);
                measures.push_back(m);
                min_ = computeMin();
            }
            double minimum() const {
                return min_;
            }
            int size() const {
                return measures.size();
            }
            void force(double value) {
                Assert(size() == maxAggregateSize);
                min_ = value;
            }
        private:
            double min_;
            std::vector<Measure> measures;
            double computeMin() const {
                std::optional<double> m;
                for(auto const & me:measures)
                {
                    if(!m || *m > me.time_by_test) {
                        m = me.time_by_test;
                    }
                }
                return *m;
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
         
         returns the count of measures added
         */
        int addMeasures(int const sz) {
            Assert(sz);
            int count = 0;
            if(0==stats.count(sz)) {
                stats[sz].include(forceTiming(sz));
                ++count;
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
                        ++count;
                        nextSz *= 2;
                    }
                }
                else {
                    stats[nextSz].include(forceTiming(nextSz));
                    ++count;
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
                        ++count;
                        prevSz /= 2;
                    }
                }
                else if(prevSz) {
                    stats[prevSz].include(forceTiming(prevSz));
                    ++count;
                }
            }
            
            return count;
        }
        
        void enforceMonotonicConstraint()
        {
            bool stable = false;
            while(!stable)
            {
                stable = true;
                backwardPass(stats,[this, &stable](auto it1, auto it2){
                    if(it1->second.minimum() > it2->second.minimum()) {
                        // assume the measure of _higher_ size is right
                        if(it1->second.size() >= maxAggregateSize) {
                            it1->second.force(it2->second.minimum());
                        }
                        else {
                            it1->second.include(forceTiming(it1->first));
                        }
                        stable = false;
                    }
                });
                forwardPass(stats,[this, &stable](auto it1, auto it2){
                    if(it1->second.minimum() > it2->second.minimum()) {
                        // assume the measure of _lower_ size is right
                        if(it2->second.size() >= maxAggregateSize) {
                            it2->second.force(it1->second.minimum());
                            return;
                        }
                        else {
                            it2->second.include(forceTiming(it2->first));
                        }
                        stable = false;
                    }
                });
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

            // to make sure there will be no hard page fault during the test,
            // we run the test (but we don't need to flush caches)
            measure_call({},
                         sz, ntests, sideEffect);
            
        start_again:

            for(int i=0; i<nConstecutiveMeasures; ++i) {
                using namespace profiling;
                
                auto microseconds = measure_call([](double & se){ polluteCache(se); },
                                                 sz, ntests, sideEffect).count();
                if(microseconds < minMicroSeconds) {
                    ntests *= 2;
                    minT.reset();
                    goto start_again;
                }
                double localTPerTest = microseconds / (1000000. * static_cast<double>(ntests));
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

        inline static CallCost cost_add_scalar_multiply { [](std::function<void(double&)> fBeforeMeasure,
                                                             int64_t sz, int64_t ntests, double & sideEffect) {
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
            if(fBeforeMeasure) {
                fBeforeMeasure(sideEffect);
            }
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
        
        inline static CallCost cost_copy { [](std::function<void(double&)> fBeforeMeasure,
                                              int64_t sz, int64_t ntests, double & sideEffect) {
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
            if(fBeforeMeasure) {
                fBeforeMeasure(sideEffect);
            }
            auto duration = measure_thread_cpu_one([&va, sz](){
                for(auto & a : va) {
                    Impl::copy(a[0].begin(),
                               a[1].begin(),
                               sz);
                }
            });
            for(auto & a : va)
            {
                using std::abs;
                sideEffect += abs(std::accumulate(a[0].begin(),
                                                  a[0].end(),
                                                  value_type{})) / static_cast<T>(a[0].size());
            }
            return duration;
        }};

        inline static CallCost cost_zero_n_raw { [](std::function<void(double&)> fBeforeMeasure,
                                                    int64_t sz, int64_t ntests, double & sideEffect) {
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
            if(fBeforeMeasure) {
                fBeforeMeasure(sideEffect);
            }
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

        inline static CallCost cost_mult_assign { [](std::function<void(double&)> fBeforeMeasure,
                                                     int64_t sz, int64_t ntests, double & sideEffect) {
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
            if(fBeforeMeasure) {
                fBeforeMeasure(sideEffect);
            }
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
        
        inline static CallCost cost_multiply { [](std::function<void(double&)> fBeforeMeasure,
                                                  int64_t sz, int64_t ntests, double & sideEffect) {
            std::vector<std::array<type, 3>> va;
            va.resize(ntests);
            for(auto & a:va) {
                int i=0;
                for(auto & v:a) {
                    ++i;
                    v.resize(2*sz);
                    setFirstReal(v, static_cast<T>(i));
                }
            }
            using namespace profiling;
            if(fBeforeMeasure) {
                fBeforeMeasure(sideEffect);
            }
            auto duration = measure_thread_cpu_one([&va, sz](){
                for(auto & a:va) {
                    Impl::multiply(a[0].data(),
                                   a[1].data(),
                                   a[2].data(),
                                   sz);
                }
            });
            for(auto & a:va) {
                sideEffect += getFirstReal(a[0]);
            }
            return duration;
        }};
        
        inline static CallCost cost_multiply_add { [](std::function<void(double&)> fBeforeMeasure,
                                                      int64_t sz, int64_t ntests, double & sideEffect) {
            std::vector<std::array<type, 3>> va;
            va.resize(ntests);
            {
                int i=0;
                for(auto & a:va) {
                    for(auto & v:a) {
                        ++i;
                        v.resize(2*sz);
                        setFirstReal(v, static_cast<T>(i));
                    }
                }
            }
            
            using namespace profiling;
            if(fBeforeMeasure) {
                fBeforeMeasure(sideEffect);
            }
            auto duration = measure_thread_cpu_one([&va, sz](){
                for(auto & a:va) {
                    Impl::multiply_add(a[0].data(),
                                       a[1].data(),
                                       a[2].data(),
                                       sz);
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

        inline static CallCost cost_fft_inverse { [](std::function<void(double&)> fBeforeMeasure,
                                                     int64_t sz, int64_t ntests, double & sideEffect) {
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
            if(fBeforeMeasure) {
                fBeforeMeasure(sideEffect);
            }
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
                    
        inline static CallCost cost_fft_forward { [](std::function<void(double&)> fBeforeMeasure,
                                                     int64_t sz, int64_t ntests, double & sideEffect) {
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
            if(fBeforeMeasure) {
                fBeforeMeasure(sideEffect);
            }
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

