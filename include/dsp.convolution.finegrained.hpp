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
        
    private:
        float cost = std::nanf("");
    };
    
    struct FinegrainedSetupParam : public Cost {
        int multiplication_group_size = 0;
    };
    
    static std::ostream& operator<<(std::ostream& os, const FinegrainedSetupParam& p)
    {
        os << "multiplication group size : " << p.multiplication_group_size;
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
            x.clear();
            x.reserve(fft_length);
            {
                y.resize(fft_length);
                it = y.begin();
            }
            
            doSetCoefficients(fft, std::move(coeffs_));
        }
        
        void setMultiplicationGroupLength(int l) {
            x.clear();
            it = y.begin();
            zero_signal(y);
            doSetMultiplicationGroupLength(l);
        }

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
            grain_number = 1;
            mult_grp_len = length;
        }
    public:
        auto getMultiplicationsGroupMaxSize() const { return mult_grp_len; }
        auto countMultiplicativeGrains() const { return 1 + (countPartitions()-1)/getMultiplicationsGroupMaxSize(); }
        auto countGrains() const { return 2 + countMultiplicativeGrains(); }
        
        int getLowestValidMultiplicationsGroupSize() const {
            // lowest valid mult_grp_len verifies:
            
            // countGrains() == getBlockSize()
            // 2 + 1 + (ffts_of_partitionned_h.size() - 1)/mult_grp_len == partition_size
            
            if(partition_size < 3) {
                // invalid configuration
                return getHighestValidMultiplicationsGroupSize();
            }
            for(int i=1;; ++i) {
                if( 2 + 1 + (ffts_of_partitionned_h.size() - 1)/i <= partition_size) {
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

            grain_number = 1;

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
        gradient_descent.setFunction( [n_frames, length_impulse, constraint, n_tests, n_channels] (int lg2_partition_size, auto & val){
            using namespace profiling;
            using namespace std;
            using namespace std::chrono;
            
            constexpr auto n_atoms_repeat = 4;
            
            if(lg2_partition_size < 0) {
                return ParamState::OutOfRange;
            }
            if(lg2_partition_size > 20) {
                throw logic_error("Gradient descent is not working?");
            }
            int const partition_size = pow2(lg2_partition_size);
            cout << "partition size : " << partition_size << endl;
            
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
                    cout << "invalid test" << endl;
                    return ParamState::OutOfRange;
                }
            }
            
            static constexpr auto index_fft = 0;
            static constexpr auto index_ifft = 1;
            array<GrainType, 2> grain_types{{ GrainType::FFT, GrainType::IFFT }};
            array<int, 2> times;
            int index = 0;
            for(auto g : grain_types)
            {
                auto prepare = [&tests, g] () { for(auto & t : tests) { t.prepare(g); } };
                auto measure = [&tests   ] () { for(auto & t : tests) { t.run();      } };
                
                // in case globals need to be initialized
                prepare(); measure();
                
                times[index] = avg(measure_n<high_resolution_clock>( n_atoms_repeat, prepare, measure));
                ++index;
            }
            
            cout << endl;
            cout << "fft  time : " << times[index_fft ] << endl;
            cout << "ifft time : " << times[index_ifft] << endl;
            
            auto max_time_fixed_grains = max(times[0], times[1]);
            float worst_grain_cost = max_time_fixed_grains;
            
            // find biggest multiplication group such that the time is just below max_fixed_grains
            
            range<int> multiplication_group_length {
                tests[0].getLowestValidMultiplicationsGroupSize(),
                tests[0].getHighestValidMultiplicationsGroupSize()
            };
            
            // in case of no spread we want to allow having multiple grains per callback.
            // But then we need to implement a gradient descent to find the sweet spot.
            // in atomic case there is no need for that because there is one less 'degree of freedom'
            // as the grain periodicity is imposed by the partition size.
            if(constraint) {
                // adjust lower bound of range to make sure that we have at most one grain every 'n_channel' audio callbacks if we use spread
                auto const min_grain_period = n_frames * n_channels;
                
                while(1) {
                    auto m = multiplication_group_length.getMin();
                    tests[0].setMultiplicationGroupLength(m);
                    if(min_grain_period <= tests[0].getGranularMinPeriod()) {
                        // Note that the caller choses the starting point of the gradient descent
                        // such that on first iteration we always break here.
                        break;
                    }
                    // grains are too close so we adjust the range
                    multiplication_group_length.setMin(multiplication_group_length.getMin() + 1);
                    if(!multiplication_group_length.empty()) {
                        continue;
                    }
                    cout << "min grain period could not be met" << endl;
                    return ParamState::OutOfRange;
                }
            }

            auto f_time_multiplications = [&tests](int mult_grp_length){
                for(auto & t : tests) {
                    t.setMultiplicationGroupLength(mult_grp_length);
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
                
                return avg(measure_n<high_resolution_clock>(n_atoms_repeat, prepare, measure));
            };
            
            auto time_mult_M = f_time_multiplications(multiplication_group_length.getMax());
            cout << "mult time for '" << multiplication_group_length.getMax() << "' (maximal) : " << time_mult_M << endl;
            if(time_mult_M <= max_time_fixed_grains) {
                // we are done because the multiplication with a maximal group length is faster
                // than fixed grains
                
                // and we prefer large periods so we chose the max group length
                multiplication_group_length.setMin(multiplication_group_length.getMax());
            }
            else {
                // the multiplication with a maximal group length is slower...
                
                auto time_mult_m = f_time_multiplications(multiplication_group_length.getMin());
                cout << "mult time for '" << multiplication_group_length.getMin() << "' (minimal) : " << time_mult_m << endl;
                if(time_mult_m >= max_time_fixed_grains) {
                    // ... and the multiplication with a minimal group is slower too
                    if(time_mult_M < time_mult_m) {
                        throw std::logic_error("more work done faster ??");
                    }
                }
                else {
                    // ... and the multiplication with a minimal group is faster,
                    // so we are in the right conditions to do a binary search

                    while(multiplication_group_length.delta() > 1) {
                        auto med = multiplication_group_length.getCenter();
                        assert(med != multiplication_group_length.getMin());
                        assert(med != multiplication_group_length.getMax());
                        
                        auto time_mult = f_time_multiplications(med);
                        cout << "mult time for '" << med << "' : " << time_mult << endl;
                        if(time_mult < max_time_fixed_grains) {
                            time_mult_m = time_mult;
                            multiplication_group_length.setMin(med);
                        }
                        else {
                            multiplication_group_length.setMax(med);
                        }
                    }
                }
                
                worst_grain_cost = time_mult_m;
                multiplication_group_length.setMax(multiplication_group_length.getMin());
            }
            
            auto optimal_group_length = multiplication_group_length.getMin();
            val.multiplication_group_size = optimal_group_length;
            
            tests[0].setMultiplicationGroupLength(optimal_group_length);
            auto period = tests[0].getGranularMinPeriod();
            
            if(constraint) {
                if(n_channels * n_frames >= period) {
                    throw std::logic_error("this should have been handled before");
                }
            }

            worst_grain_cost /= n_tests;
            // worst_grain_cost == 'one computation'
            
            auto n_max_computes_per_callback = n_frames / period;
            if(n_frames != n_max_computes_per_callback * period) {
                // in the worst case, we have one more
                ++ n_max_computes_per_callback;
            }
            if(constraint) {
                if(n_max_computes_per_callback != 1) {
                    throw logic_error("the constraint ensures that the number of"
                                      " computes per callback is 1/n_channels on average");
                }
                // n_frames is small enough and partition_size is big enough so that
                // there is enough "room" to spread the computes of different channels over different callback calls,
                // provided we "phase" the different partitionned convolutions correctly.
                // Hence we take this advantage into account here:
                worst_grain_cost /= n_channels;
            }
            
            worst_grain_cost *= n_max_computes_per_callback;
            // worst_grain_cost == 'worst computation time over one callback'
            
            worst_grain_cost /= n_frames;
            // worst_grain_cost == 'worst computation time over one callback, averaged per frame'
            
            val.setCost(worst_grain_cost);
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
        
        return gradient_descent.findMinimum(n_iterations,
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
        /* timings have random noise, so iterating helps having a better precision */
        constexpr auto n_iterations = 30;
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
