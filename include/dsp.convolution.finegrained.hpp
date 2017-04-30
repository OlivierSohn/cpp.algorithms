/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule
{
    enum class GrainType {
        FFT,
        IFFT,
        MultiplicationGroup
    };
    
    struct Cost {
        void setCost(float c) { cost = c; }
        
        operator float() const { return cost; }
        operator float&() { return cost; }
        
        float getCost() const { return cost; }
    private:
        float cost = std::numeric_limits<float>::quiet_NaN();
    };
    
    struct GrainsCosts {
        float fft, ifft, mult;
    };
    
    static std::ostream& operator<<(std::ostream& os, const GrainsCosts& g)
    {
        using namespace std;
        os << "grain 'fft'  : " << g.fft << endl;
        os << "grain 'ifft' : " << g.ifft << endl;
        os << "grain 'mult' : " << g.mult << endl;
        return os;
    }

    struct FinegrainedSetupParam : public Cost {
        void setPhase(int ph) { phase = ph; }
        void setGrainsCosts(GrainsCosts gcosts) { grains_costs = gcosts; }
        
        int multiplication_group_size = 0;
        int phase = 0;
        GrainsCosts grains_costs;
    };
    
    static std::ostream& operator<<(std::ostream& os, const FinegrainedSetupParam& p)
    {
        using namespace std;
        os
        << "phase : " << p.phase << endl
        << p.grains_costs
        << "multiplication group size : " << p.multiplication_group_size << endl;
        return os;
    }
    
    template <typename Parent>
    struct FinegrainedFFTConvolutionBase : public Parent {
        using T = typename Parent::FPT;
        using FPT = T;
        using Tag = typename Parent::FFTTag;
        
        
        using SetupParam = FinegrainedSetupParam;
        
        void applySetup(SetupParam const & p) {
            setMultiplicationGroupLength(p.multiplication_group_size);
        }
        
        static constexpr auto zero_signal = fft::RealSignal_<Tag, FPT>::zero;
        static constexpr auto copy = fft::RealSignal_<Tag, FPT>::copy;
        static constexpr auto get_signal = fft::RealSignal_<Tag, FPT>::get_signal;
        static constexpr auto add_scalar_multiply = fft::RealSignal_<Tag, FPT>::add_scalar_multiply;
        using RealSignal = typename fft::RealSignal_<Tag, FPT>::type;
        using Algo = typename fft::Algo_<Tag, FPT>;
        using Contexts = fft::Contexts_<Tag, FPT>;
        
        using Parent::get_fft_length;
        using Parent::getBlockSize;
        using Parent::doSetCoefficients;
        
        using Parent::getGranularMinPeriod;
        using Parent::doSetMultiplicationGroupLength;
        using Parent::isValid;
        using Parent::countPartitions;
        using Parent::countGrains;
        using Parent::getGrainNumber;
        using Parent::increment_grain;
        using Parent::compute_x_fft;
        using Parent::do_some_multiply_add;
        using Parent::get_multiply_add_result;
        
        static constexpr bool is_atomic = false;
        
        FinegrainedFFTConvolutionBase() {
            it = y.end();
        }
        
        void setCoefficients(a64::vector<T> coeffs_) {
            
            if(coeffs_.size() < 2) {
                coeffs_.resize(2); // avoid ill-formed cases
            }
            auto const N = coeffs_.size();
            auto const fft_length = get_fft_length(N);
            fft.setContext(Contexts::getInstance().getBySize(fft_length));
            
            result.resize(fft_length);
            x.reserve(fft_length);
            {
                y.resize(fft_length);
            }
            
            doSetCoefficients(fft, std::move(coeffs_));
            reset_states();
        }
        
        void setMultiplicationGroupLength(int l) {
            doSetMultiplicationGroupLength(l);
            reset_states();
        }
        
    private:
        void reset_states() {
            Parent::reset_base_states();
            x.clear();
            zero_signal(y);
            it = y.begin();
            grain_counter = 0;
        }
        
    public:
        void fastForwardToComputation(GrainType t, T val = 1) {
            switch(t) {
                case GrainType::FFT:
                    while(x.size() != getBlockSize()-1) {
                        step(val);
                    }
                    break;
                case GrainType::IFFT:
                    while(getGrainNumber() != countGrains() - 1) {
                        step(val);
                    }
                    while(grain_counter != getGranularMinPeriod() - 1) {
                        step(val);
                    }
                    break;
                case GrainType::MultiplicationGroup:
                    while(!x.empty()) {
                        step(val);
                    }
                    while(grain_counter != getGranularMinPeriod() - 1) {
                        step(val);
                    }
                    break;
            }
            assert(willComputeNextStep());
        }
        
        bool willComputeNextStep() const {
            return (x.size() == getBlockSize()-1) || (grain_counter+1 == getBlockSize()/countGrains());
        }
        
        void step(T val) {
            assert(isValid());
            x.emplace_back(val);
            
            auto const block_size = getBlockSize();
            auto const sz = x.size();
            if(unlikely(sz == block_size)) {
                // pad x
                
                x.resize(get_fft_length());
                
                compute_x_fft(fft, x);
                
                x.clear();
                
                // in the same "grain of computation" we do the following.
                // if it is too costly we must delay the computation
                // in another grain and have a second y vector
                
                auto factor = 1 / (Algo::scale * Algo::scale * static_cast<T>(get_fft_length()));
                
                auto it_res = result.begin();
                auto it_y = y.begin();
                auto it_y_prev = it_y + block_size;
                
                // y = mix first part of result with second part of previous result
                //
                // 'first part of y' = factor * ('second part of y' + 'first part of result')
                add_scalar_multiply(it_y, /* = */
                                    /* ( */ it_res, /* + */ it_y_prev /* ) */, /* x */ factor,
                                    block_size);
                
                // store second part of result for later
                //
                // 'second part of y' = 'second part of result'
                copy(it_y   + block_size,
                     it_res + block_size,
                     block_size);
                
                // reset 'it' so that the results are accessible in get() method
                assert(it == y.begin() + block_size-1); // make sure 'rythm is good', i.e we exhausted the first half of the y vector
                it = y.begin();
                grain_counter = 0;
                increment_grain();
                return;
            }
            ++grain_counter;
            
            auto const n_grains = countGrains();
            auto const granularity = block_size/n_grains;
            assert(granularity >= grain_counter);
            if(grain_counter == granularity) {
                grain_counter = 0;
                auto cur_grain = getGrainNumber();
                assert(cur_grain <= n_grains);
                if( cur_grain < n_grains - 1 ) {
                    do_some_multiply_add();
                }
                else if(cur_grain == n_grains - 1) {
                    fft.inverse(get_multiply_add_result(),
                                result,
                                get_fft_length());
                    increment_grain();
                }
                else {
                    // spread is not optimal
                }
            }
            ++it;
        }
        
        T get() {
            assert(it < y.begin() + getBlockSize());
            assert(it >= y.begin());
            return get_signal(*it);
        }
        
    private:
        int grain_counter = 0;
        Algo fft;
        RealSignal x, y, result;
        typename decltype(y)::iterator it;
    };
    
    template <typename T, typename Tag>
    struct FinegrainedPartitionnedFFTConvolutionCRTP {
        using FPT = T;
        using FFTTag = Tag;
        
        using RealSignal = typename fft::RealSignal_<Tag, FPT>::type;
        using Signal_value_type = typename fft::RealSignal_<Tag, FPT>::value_type;
        
        using CplxFreqs = typename fft::RealFBins_<Tag, FPT>::type;
        static constexpr auto zero = fft::RealFBins_<Tag, FPT>::zero;
        static constexpr auto multiply_add = fft::RealFBins_<Tag, FPT>::multiply_add;
        
        using Algo = typename fft::Algo_<Tag, FPT>;
        
        auto get_fft_length() const { assert(partition_size > 0); return 2 * partition_size; }
        auto get_fft_length(int) const { return get_fft_length(); }
        
        bool empty() const { return ffts_of_partitionned_h.empty(); }
        
        auto getBlockSize() const { return partition_size; }
        auto getLatency() const { return 2*partition_size; }
        auto getGranularMinPeriod() const { return getBlockSize() / countGrains(); }
        bool isValid() const { return countGrains() <= getBlockSize(); }
        
        int countPartitions() const { return ffts_of_partitionned_h.size(); }
        
    protected:
        void doSetMultiplicationGroupLength(int length) {
            mult_grp_len = length;
        }
        
        void reset_base_states() {
            grain_number = 1;
        }
    public:
        auto getMultiplicationsGroupMaxSize() const { return mult_grp_len; }
        auto countMultiplicativeGrains() const { return 1 + (countPartitions()-1)/getMultiplicationsGroupMaxSize(); }
        static constexpr auto countNonMultiplicativeGrains() { return 2; }
        auto countGrains() const { return countNonMultiplicativeGrains() + countMultiplicativeGrains(); }
        
        int getLowestValidMultiplicationsGroupSize() const {
            // lowest valid mult_grp_len verifies:
            
            // countGrains() == getBlockSize()
            // 2 + 1 + (ffts_of_partitionned_h.size() - 1)/mult_grp_len == partition_size
            
            auto constexpr min_number_grains = countNonMultiplicativeGrains() + 1;
            auto diff = partition_size - min_number_grains;
            if(diff < 0) {
                // invalid configuration
                return getHighestValidMultiplicationsGroupSize();
            }
            assert(diff >=0 );
            for(int i=1;; ++i) {
                if( (ffts_of_partitionned_h.size() - 1)/i <= diff) {
                    return i;
                }
            }
        }
        
        int getHighestValidMultiplicationsGroupSize() const { return countPartitions(); }
        
        void set_partition_size(int sz) {
            assert(sz > 0);
            partition_size = sz;
            assert(is_power_of_two(sz));
        }
        
        void doSetCoefficients(Algo const & fft, a64::vector<T> coeffs_) {
            
            auto const n_partitions = [&coeffs_, partition_size = this->partition_size](){
                auto const N = coeffs_.size();
                auto n_partitions = N/partition_size;
                if(n_partitions * partition_size != N) {
                    // one partition is partial...
                    assert(n_partitions * partition_size < N);
                    ++n_partitions;
                    // ... pad it with zeros
                    coeffs_.resize(n_partitions * partition_size, {});
                }
                return n_partitions;
            }();
            
            ffts_of_delayed_x.reset();
            ffts_of_delayed_x.resize(n_partitions);
            ffts_of_partitionned_h.resize(n_partitions);
            
            auto const fft_length = get_fft_length();
            
            for(auto & fft_of_partitionned_h : ffts_of_partitionned_h) {
                fft_of_partitionned_h.resize(fft_length);
            }
            for(auto & fft_of_delayed_x : ffts_of_delayed_x) {
                fft_of_delayed_x.resize(fft_length);
            }
            
            work.resize(fft_length);
            
            // compute fft of padded impulse response
            
            auto it_coeffs = coeffs_.begin();
            {
                RealSignal coeffs_slice(fft_length, Signal_value_type(0)); // initialize with zeros (second half is padding)
                for(auto & fft_of_partitionned_h : ffts_of_partitionned_h) {
                    auto end_coeffs = it_coeffs + partition_size;
                    assert(end_coeffs <= coeffs_.end());
                    auto slice_it = coeffs_slice.begin();
                    for(;it_coeffs != end_coeffs; ++it_coeffs, ++slice_it) {
                        using RealT = typename RealSignal::value_type;
                        *slice_it = RealT(*it_coeffs);
                    }
                    
                    // coeffs_slice is padded with 0, because it is bigger than partition_size
                    // and initialized with zeros.
                    fft.forward(coeffs_slice, fft_of_partitionned_h, fft_length);
                }
            }
            assert(it_coeffs == coeffs_.end());
        }
        
    protected:
        
        int getGrainNumber() const { return grain_number; }
        
        void compute_x_fft(Algo const & fft, RealSignal const & x) {
            assert(grain_number == countGrains());
            grain_number = 0;
            auto & oldest_fft_of_delayed_x = *ffts_of_delayed_x.cycleEnd();
            auto const fft_length = get_fft_length();
            assert(fft_length == oldest_fft_of_delayed_x.size());
            fft.forward(x, oldest_fft_of_delayed_x, fft_length);
            ffts_of_delayed_x.advance();
        }
        
        void do_some_multiply_add() {
            auto const M = getMultiplicationsGroupMaxSize();
            auto const offset_base = M * (grain_number - 1);
            assert(offset_base >= 0);
            if(0 == offset_base) {
                zero(work);
            }
            
            assert(offset_base < ffts_of_partitionned_h.size());
            for(auto i=0; i < M; ++i) {
                auto offset = offset_base + i;
                if(offset == ffts_of_partitionned_h.size()) {
                    break;
                }
                auto it_fft_of_partitionned_h = ffts_of_partitionned_h.begin() + offset;
                assert(it_fft_of_partitionned_h < ffts_of_partitionned_h.end());
                
                auto const & fft_of_delayed_x = ffts_of_delayed_x.get_backward(offset);
                
                multiply_add(work                /*   +=   */,
                             fft_of_delayed_x,   /*   x   */   *it_fft_of_partitionned_h);
            }
            increment_grain();
        }
        
        auto const & get_multiply_add_result() const {
            return work;
        }
        
        void increment_grain() {
            ++grain_number;
        }
        
    private:
        int mult_grp_len = 0;
        int partition_size = -1;
        int grain_number = 0;
        cyclic<CplxFreqs> ffts_of_delayed_x;
        std::vector<CplxFreqs> ffts_of_partitionned_h;
        
        CplxFreqs work;
    };
    
    template <typename T, typename FFTTag = fft::Fastest>
    using FinegrainedPartitionnedFFTConvolution = FinegrainedFFTConvolutionBase< FinegrainedPartitionnedFFTConvolutionCRTP<T, FFTTag> >;
    
    /*
     input parameters :
     - n_frames_audio_cb, n_channels, with_spread (those 3 can be reduced to 'equivalent_n_frames_cb')
     - impulse response length
     
     output parameters:
     - lg(partition size)
     
     1D - Gradient descent according to cost 'max grain computation time' with variable parameters 'lg(partition_size)'
     deducing 'number of multiplications per grain' by finding the parameter that leads to grain computations time just below max(fft, ifft),
     with the constraint that one computation at most occurs per 'equivalent' audio callback.
     */
    template<typename NonAtomicConvolution, typename SetupParam = typename NonAtomicConvolution::SetupParam, typename GradientDescent = GradientDescent<SetupParam>>
    int get_lg2_optimal_partition_size_for_nonatomic_convolution(GradientDescent & gradient_descent,
                                                                 int n_iterations,
                                                                 int n_channels,
                                                                 int n_frames,
                                                                 int length_impulse,
                                                                 bool constraint,
                                                                 SetupParam & min_val,
                                                                 int n_tests) {
        //std::cout << "main thread: " << std::endl;
        //thread::logSchedParams();
        
        gradient_descent.setFunction( [n_frames, length_impulse, constraint, n_tests, n_iterations, n_channels] (int lg2_partition_size, auto & val){
            using namespace profiling;
            using namespace std;
            using namespace std::chrono;
            
            constexpr auto n_atoms_repeat = 1;
            
            if(lg2_partition_size < 0) {
                return ParamState::OutOfRange;
            }
            if(lg2_partition_size > 20) {
                throw logic_error("Gradient descent is not working?");
            }
            int const partition_size = pow2(lg2_partition_size);
            //            cout << "partition size : " << partition_size << endl;
            
            struct Test {
                using T = typename NonAtomicConvolution::FPT;
                Test(size_t partition_size, int length_impulse) {
                    pfftcv.set_partition_size(partition_size);
                    pfftcv.setCoefficients(a64::vector<T>(length_impulse));
                    
                    // the value is not very important it will be overriden later on,
                    // but we need to set something valid
                    pfftcv.setMultiplicationGroupLength(pfftcv.getHighestValidMultiplicationsGroupSize());
                }
                
                bool isValid() const { return pfftcv.isValid(); }
                
                int getGranularMinPeriod() const { return pfftcv.getGranularMinPeriod(); }
                int countMultiplicativeGrains() const { return pfftcv.countMultiplicativeGrains(); }
                int getHighestValidMultiplicationsGroupSize() const { return pfftcv.getHighestValidMultiplicationsGroupSize(); }
                int getLowestValidMultiplicationsGroupSize() const { return pfftcv.getLowestValidMultiplicationsGroupSize(); }
                
                void setMultiplicationGroupLength(int mult_grp_length) {
                    pfftcv.setMultiplicationGroupLength(mult_grp_length);
                }
                
                void prepare(GrainType g) {
                    pfftcv.fastForwardToComputation(g);
                }
                void run() {
                    assert(pfftcv.willComputeNextStep());
                    pfftcv.step(1.f);
                    pfftcv.get();
                }
            private:
                NonAtomicConvolution pfftcv;
            };
            
            // prepare tests
            
            vector<Test> tests;
            tests.reserve(n_tests);
            for(int i=0; i<n_tests;++i) {
                tests.emplace_back(partition_size, length_impulse);
                if(!tests[i].isValid()) {
                    //cout << "invalid test" << endl;
                    return ParamState::OutOfRange;
                }
            }
            
            constexpr auto n_non_multiplicative_grains = NonAtomicConvolution::countNonMultiplicativeGrains();
            static_assert(2 == n_non_multiplicative_grains, "");
            static constexpr auto index_fft = 0;
            static constexpr auto index_ifft = 1;
            array<GrainType, n_non_multiplicative_grains> grain_types{{ GrainType::FFT, GrainType::IFFT }};
            array<float, n_non_multiplicative_grains> times;
            int index = 0;
            for(auto g : grain_types)
            {
                auto prepare = [&tests, g] () { for(auto & t : tests) { t.prepare(g); } };
                auto measure = [&tests   ] () { for(auto & t : tests) { t.run();      } };
                
                // in case globals need to be initialized
                prepare(); measure();
                
                times[index] = min_(measure_n<high_resolution_clock>(n_atoms_repeat,
                                                                     prepare,
                                                                     measure)) / static_cast<float>(tests.size());
                ++index;
            }
            
            //            cout << endl;
            //            cout << "fft  time : " << times[index_fft ] << endl;
            //            cout << "ifft time : " << times[index_ifft] << endl;
            
            struct PhasedCost : public Cost {
                int phase = 0; // in frames
                GrainsCosts grains_costs;
            };

            struct CostEvaluator {
                array<float, n_non_multiplicative_grains> fft_times;
                int n_audio_cb_frames;
                int n_channels;
                bool constraint;
                
                void evaluate(float multiplication_grain_time, int n_multiplicative_grains, int period,
                               PhasedCost & result) const {
                    result.phase = 0;
                    
                    auto max_n_grains_per_cb = n_audio_cb_frames / period;
                    if(max_n_grains_per_cb * period != n_audio_cb_frames) {
                        ++max_n_grains_per_cb; // worst case, not avg
                    }
                    
                    cyclic<float> grains_costs(n_multiplicative_grains + n_non_multiplicative_grains,
                                               multiplication_grain_time); // initialize with multiplicative grain times
                    for(auto t : fft_times) {
                        grains_costs.feed(t);
                    }
                    
                    float cost = computeMaxSlidingSum(grains_costs,
                                                      max_n_grains_per_cb);
                    
                    if(constraint) {
                        assert(n_channels >= 2);
                        auto n_min_empty_cb_between_consecutive_grains = -1 + period / n_audio_cb_frames;
                        if(n_min_empty_cb_between_consecutive_grains >= n_channels - 1) {
                            // easy case : there is enough room between grains to evenly distribute all channels
                            result.phase = period / n_channels;
                        }
                        else {
                            // harder case: we need to go more in detail, and find the phase that minimizes
                            // the worst callback cost.

                            cost *= n_channels;
                            // now cost is the 'phase == 0' cost
                            result.phase = 0;
                            
                            cyclic<float> phased_grains_costs;
                            for(int phase = 1; phase < grains_costs.size(); ++phase) {
                                compute_phased_sum(grains_costs,
                                                   phase,
                                                   n_channels,
                                                   phased_grains_costs);
                                
                                auto phased_cost = computeMaxSlidingSum(phased_grains_costs,
                                                                        max_n_grains_per_cb);
                                if(phased_cost < cost) {
                                    cost = phased_cost;
                                    result.phase = phase;
                                }
                            }
                            
                            // convert phase units from "grain periods" to "frames"
                            result.phase *= period;
                        }
                        // cost is now the cost of each channel
                        // but cost should be per sample, not per frame, so
                        // we divide by the number of channels
                        cost /= static_cast<float>(n_channels);
                    }
                    
                    cost /= n_audio_cb_frames;
                    // cost == 'worst computation time over one callback, averaged per sample'
                    
                    result.grains_costs.fft  = fft_times[index_fft];
                    result.grains_costs.ifft = fft_times[index_ifft];
                    result.grains_costs.mult = multiplication_grain_time;

                    result.setCost(cost);
                }
            } cost_evaluator{times, n_frames, n_channels, constraint};
            
            RangedGradientDescent<PhasedCost> rgd([ &cost_evaluator, &tests ](int multiplication_group_size, auto & cost) {
                // compute multiplication time for the group
                
                for(auto & t : tests) {
                    t.setMultiplicationGroupLength(multiplication_group_size);
                }
                
                if(!tests[0].isValid()) {
                    return ParamState::OutOfRange;
                }
                
                auto prepare = [&tests]() {
                    for(auto & t : tests) {
                        t.prepare(GrainType::MultiplicationGroup);
                    }
                };
                auto measure = [&tests](){
                    for(auto & t : tests) {
                        t.run();
                    }
                };
                
                auto multiplication_grain_time = min_(measure_n<high_resolution_clock>(n_atoms_repeat,
                                                                                       prepare,
                                                                                       measure)) / static_cast<float>(tests.size());
                
                cost_evaluator.evaluate(multiplication_grain_time,
                                        tests[0].countMultiplicativeGrains(),
                                        tests[0].getGranularMinPeriod(),
                                        cost);
                
                //cout
                //<< "mult time for group size '" << multiplication_group_size << "' : " << multiplication_grain_time
                //<< " cost : '" << cost << "'" << endl;
                
                return ParamState::Ok;
            });
            
            range<int> const multiplication_group_length {
                tests[0].getLowestValidMultiplicationsGroupSize(),
                tests[0].getHighestValidMultiplicationsGroupSize()
            };
            
            PhasedCost phased_cost;
            val.multiplication_group_size = rgd.findLocalMinimum(n_iterations, multiplication_group_length, phased_cost);
            val.setCost(phased_cost.getCost());
            val.setGrainsCosts(phased_cost.grains_costs);
            val.setPhase(phased_cost.phase);
            
            constexpr auto debug = false;
            if(debug) {
                rgd.plot();
                rgd.make_exhaustive(multiplication_group_length);
                rgd.plot();
            }
            
            //            cout
            //            << "optimal group size : " << val.multiplication_group_size
            //            << " cost : '" << cost << "'" << endl;
            
            return ParamState::Ok;
        });
        
        auto start_lg2_partition = 7;
        // to ensure that the constraint is met in first try
        if(constraint) {
            start_lg2_partition =
            2 // because we have 3 or 4 grains at least (2^2)
            + 1 + power_of_two_exponent(n_channels * n_frames); // criteria for atomic version of convolution
        }
        else {
            start_lg2_partition =
            2 // because we have 3 or 4 grains at least (2^2)
            + 1 + power_of_two_exponent(n_frames); // criteria for atomic version of convolution
        }
        
        return gradient_descent.findLocalMinimum(n_iterations,
                                                 start_lg2_partition,
                                                 min_val);
    }
    
    template<typename NonAtomicConvolution, typename SetupParam = typename NonAtomicConvolution::SetupParam, typename GradientDescent = GradientDescent<SetupParam>>
    int get_optimal_partition_size_for_nonatomic_convolution(GradientDescent & gd,
                                                             int n_channels,
                                                             bool with_spread,
                                                             int n_audiocb_frames,
                                                             int length_impulse,
                                                             SetupParam & value )
    {
        constexpr auto n_iterations = 1;
        constexpr auto n_tests = 1;
        int lg2_part_size = get_lg2_optimal_partition_size_for_nonatomic_convolution<NonAtomicConvolution>(gd,
                                                                                                           n_iterations,
                                                                                                           n_channels,
                                                                                                           n_audiocb_frames,
                                                                                                           length_impulse,
                                                                                                           with_spread,
                                                                                                           value,
                                                                                                           n_tests);
        
        return pow2(lg2_part_size);
    }
    
    template<typename T>
    struct PartitionAlgo< FinegrainedPartitionnedFFTConvolution<T> > {
        using NonAtomicConvolution = FinegrainedPartitionnedFFTConvolution<T>;
        using SetupParam = typename NonAtomicConvolution::SetupParam;
        using PartitionningSpec = PartitionningSpec<SetupParam>;
        using PartitionningSpecs = PartitionningSpecs<SetupParam>;
        
        static PartitionningSpecs run(int n_channels, int n_audio_cb_frames, int size_impulse_response) {
            assert(n_channels > 0);
            PartitionningSpecs res;
            {
                auto & spec = res.without_spread;
                spec.size = get_optimal_partition_size_for_nonatomic_convolution<NonAtomicConvolution>(spec.gd,
                                                                                                       n_channels,
                                                                                                       false,
                                                                                                       n_audio_cb_frames,
                                                                                                       size_impulse_response,
                                                                                                       spec.cost );
            }
            
            if(n_channels > 1) {
                auto & spec = res.with_spread;
                spec.size = get_optimal_partition_size_for_nonatomic_convolution<NonAtomicConvolution>(spec.gd,
                                                                                                       n_channels,
                                                                                                       true,
                                                                                                       n_audio_cb_frames,
                                                                                                       size_impulse_response,
                                                                                                       spec.cost );
            }
            
            return std::move(res);
        }
    };
}
