/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule
{    
    /*
     * cf. https://en.wikipedia.org/wiki/Overlap%E2%80%93add_method#math_Eq.1
     */
    template <typename Parent>
    struct FFTConvolutionBase : public Parent {
        using T = typename Parent::FPT;
        using FPT = T;
        using Tag = typename Parent::FFTTag;
        static constexpr auto copy = fft::RealSignal_<Tag, FPT>::copy;
        static constexpr auto get_signal = fft::RealSignal_<Tag, FPT>::get_signal;
        static constexpr auto add_scalar_multiply = fft::RealSignal_<Tag, FPT>::add_scalar_multiply;
        using RealSignal = typename fft::RealSignal_<Tag, FPT>::type;
        using Algo = typename fft::Algo_<Tag, FPT>;
        using Contexts = fft::Contexts_<Tag, FPT>;
        
        using Parent::get_fft_length;
        using Parent::getComputationPeriodicity;
        using Parent::doSetCoefficients;
        using Parent::compute_convolution;
        
        FFTConvolutionBase() {
            it = y.end();
        }
        
        void setCoefficients(std::vector<T> coeffs_) {
            
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
                it = y.begin();
            }

            doSetCoefficients(fft, std::move(coeffs_));
        }
        
        void step(T val) {
            x.emplace_back(val);
            
            auto N = getComputationPeriodicity();
            if(unlikely(x.size() == N)) {
                // pad x
                
                x.resize(get_fft_length());
                
                auto const & frequencies = compute_convolution(fft, x);

                fft.inverse(frequencies, result, get_fft_length());

                auto factor = 1 / (Algo::scale * Algo::scale * static_cast<T>(get_fft_length()));
                
                auto it_res = result.begin();
                auto it_y = y.begin();
                auto it_y_prev = it_y + N;
                
                // y = mix first part of result with second part of previous result
                //
                // 'first part of y' = factor * ('second part of y' + 'first part of result')
                add_scalar_multiply(it_y, /* = */
                                    /* ( */ it_res, /* + */ it_y_prev /* ) */, /* x */ factor,
                                    N);
                
                // store second part of result for later
                //
                // 'second part of y' = 'second part of result'
                copy(it_y   + N,
                     it_res + N,
                     N);
                
                // reset 'it' so that the results are accessible in get() method
                assert(it == y.begin() + N-1); // make sure 'rythm is good', i.e we exhausted the first half of the y vector
                it = y.begin();
                
                // clear 'x'
                
                x.clear();
            }
            else {
                ++it;
            }
        }
        
        T get() {
            assert(it < y.begin() + getComputationPeriodicity());
            assert(it >= y.begin());
            return get_signal(*it);
        }
        
    private:
        Algo fft;
        RealSignal x, y, result;
        typename decltype(y)::iterator it;
    };
    
    template <typename T, typename Tag>
    struct FFTConvolutionCRTP {
        using FPT = T;
        using FFTTag = Tag;
        
        using RealSignal = typename fft::RealSignal_<Tag, FPT>::type;
        static constexpr auto makeRealSignal = fft::RealSignal_<Tag, FPT>::make;

        using CplxFreqs = typename fft::RealFBins_<Tag, FPT>::type;
        static constexpr auto mult_assign = fft::RealFBins_<Tag, FPT>::mult_assign;
        
        using Algo = typename fft::Algo_<Tag, FPT>;

        auto getComputationPeriodicity() const { return N; }
        
        auto countPartitions() const { return 1; }
        
        bool empty() const { return fft_of_h.empty(); }
        
        auto get_fft_length(int N) {
            auto N_nonzero_y = 2 * N - 1;
            return ceil_power_of_two(N_nonzero_y);
        }
        
        auto get_fft_length() {
            return fft_of_x.size();
        }
        
        void doSetCoefficients(Algo const & fft, std::vector<T> coeffs_) {
            
            N = coeffs_.size();
            
            auto fft_length = get_fft_length(N);
            
            fft_of_h.resize(fft_length);
            fft_of_x.resize(fft_length);
            
            // pad impulse response with 0
            
            coeffs_.resize(fft_length, {});
            
            // compute fft of padded impulse response
            auto coeffs = makeRealSignal(std::move(coeffs_));
            fft.forward(coeffs, fft_of_h, fft_length);
        }
        
    protected:

        auto const & compute_convolution(Algo const & fft, RealSignal const & x) {
            // do fft of x
            
            fft.forward(x, fft_of_x, get_fft_length());
            
            mult_assign(fft_of_x, fft_of_h);            
            
            return fft_of_x;
        }
      
    private:
        int N;
        
        CplxFreqs fft_of_h, fft_of_x;
    };
    
    /*
     * Partitionned convolution, cf. http://www.ericbattenberg.com/school/partconvDAFx2011.pdf
     *
     * The impulse response h is split in parts of equal length h1, h2, ... hn
     *
     *                     FFT(h1)
     *                        |
     *       +-----+    +-----v-----+  +-----+  +------+
     * x +---> FFT +-+--> cplx mult +--> Add +--> IFFT +--> y
     *       +-----+ |  +-----------+  +-^---+  +------+
     *          +----v---+               |
     *          |Delay(N)| FFT(h2)       |
     *          +----+---+    |          |
     *               |  +-----v-----+    |
     *               +->+ cplx mult +----+
     *               |  +-----------+    |
     *               .         .         .
     *               .         .         .
     *               .                   .
     *               |                   |
     *          +----v---+               |
     *          |Delay(N)| FFT(hn)       |
     *          +----+---+    |          |
     *               |  +-----v-----+    |
     *               +->+ cplx mult +----+
     *                  +-----------+
     */
    
    
    int get_lg2_optimal_partition_size(GradientDescent & gd,
                                       int n_iterations,
                                       int n_channels,
                                       int n_frames,
                                       int length_impulse,
                                       bool constraint,
                                       float & min_val,
                                       int n_tests);
    
    static inline int get_optimal_partition_size(GradientDescent & gd,
                                                 int n_channels,
                                                 bool with_spread,
                                                 int n_audiocb_frames,
                                                 int length_impulse,
                                                 float & value )
    {
        /* timings have random noise, so iterating helps having a better precision */
        constexpr auto n_iterations = 30;
        constexpr auto n_tests = 1;
        int lg2_part_size = get_lg2_optimal_partition_size(gd,
                                                           n_iterations,
                                                           n_channels,
                                                           n_audiocb_frames,
                                                           length_impulse,
                                                           with_spread,
                                                           value,
                                                           n_tests);
        
        return pow2(lg2_part_size);
    }
    
    enum class Autotune {
        Yes,
        No
    };

    template <typename T, typename Tag>
    struct PartitionnedFFTConvolutionCRTP {
        using FPT = T;
        using FFTTag = Tag;

        using RealSignal = typename fft::RealSignal_<Tag, FPT>::type;
        
        using CplxFreqs = typename fft::RealFBins_<Tag, FPT>::type;
        static constexpr auto zero = fft::RealFBins_<Tag, FPT>::zero;
        static constexpr auto multiply_add = fft::RealFBins_<Tag, FPT>::multiply_add;
        
        using Algo = typename fft::Algo_<Tag, FPT>;

        auto get_fft_length() { assert(partition_size > 0); return 2 * partition_size; }
        auto get_fft_length(int) { return get_fft_length(); }
        
        bool empty() const { return ffts_of_partitionned_h.empty(); }

        auto getComputationPeriodicity() { return partition_size; }
        
        auto countPartitions() const { return ffts_of_partitionned_h.size(); }
        
        void set_partition_size(int sz) {
            assert(sz > 0);
            partition_size = sz;
            assert(is_power_of_two(sz));
        }

        void doSetCoefficients(Algo const & fft, std::vector<T> coeffs_) {
            
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
                RealSignal coeffs_slice(fft_length, {}); // initialize with zeros (second half is padding)
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

        auto const & compute_convolution(Algo const & fft, RealSignal const & x)
        {
            // do fft of x
            {
                auto & oldest_fft_of_delayed_x = *ffts_of_delayed_x.cycleEnd();
                auto const fft_length = get_fft_length();
                assert(fft_length == oldest_fft_of_delayed_x.size());
                fft.forward(x, oldest_fft_of_delayed_x, fft_length);
                ffts_of_delayed_x.advance();
            }
            
            auto it_fft_of_partitionned_h = ffts_of_partitionned_h.begin();
            
            zero(work);
            
            ffts_of_delayed_x.for_each_bkwd( [this, &it_fft_of_partitionned_h] (auto const & fft_of_delayed_x) {
                assert(it_fft_of_partitionned_h < ffts_of_partitionned_h.end());
                
                // work += fft_of_delayed_x * fft_of_partitionned_h
                multiply_add(work, fft_of_delayed_x, *it_fft_of_partitionned_h);
                
                ++ it_fft_of_partitionned_h;
            });
            
            return work;
        }

    private:
        int partition_size = -1;
        cyclic<CplxFreqs> ffts_of_delayed_x;
        std::vector<CplxFreqs> ffts_of_partitionned_h;
        
        CplxFreqs work;
    };
    
    /*
    
     Notations for complexity:
     
     H : length of impulse response
     
     A : number of frames computed during one audio callback
     
     PART_N : Number of partitions
     PART_L : Length of one partition
     ( H == PART_N * PART_L )
     
     The computations are based on the fact that an fft of an S-long signal costs S*lg(S)
     */

    // for convolution reverbs, most of the time we have A << H
    /*
     * runtime complexity:
     *
     *   every H frames  ............... O( H * lg(H) )
     *
     *   every frame  .................. O( lg(H) )       [amortized]
     *
     *   when 'H < A':
     *     worst audio callback call ... O( A * lg(H) )    (= A/H * O( H * lg(H) ))
     *
     *   when 'A < H ':
     *     worst audio callback call ... O( H * lg(H) )
     *
     * optimization : H and A are fixed so we cannot optimize this algorithm
     */
    template <typename T, typename FFTTag = fft::Fastest>
    using FFTConvolution = FFTConvolutionBase< FFTConvolutionCRTP<T, FFTTag> >;

    /*
     * runtime complexity:
     *
     *   every PART_L frames                :   O( PART_L * (PART_N + lg(PART_L) ) )
     *
     *   every frame                        :   O( PART_N + lg(PART_L) )      [amortized]
     *
     *   with 'PART_L < A':
     *     worst audio callback call ... O(     A  * (part_N + lg(PART_L)) )       (= A/PART_L *  O( PART_L * (PART_N + lg(PART_L) ) )
     *
     *   with 'A < PART_L ' 
     *     worst audio callback call ... O( PART_L * (PART_N + lg(PART_L) ) )
     *                                 = O( PART_L * (H/PART_L + lg(PART_L) ) )
     *                                 = O( H + PART_L * lg(PART_L) ) )
     *
     * optimization : PART_L is not fixed so we can optimize this algorithm by trying different powers of 2
     *                also we can optimize more globally, taking into account that we have one reverb per channel:
     *                when N_Channel * A < PART_L we can distribute the computes over the different callbacks calls, provided the "phases"
     *                of the different algorithms are well-spaced.
     *
     */
    template <typename T, typename FFTTag = fft::Fastest>
    using PartitionnedFFTConvolution = FFTConvolutionBase< PartitionnedFFTConvolutionCRTP<T, FFTTag> >;
    
    
    struct PartitionningSpec {
        int size =  -1;
        float avg_time_per_sample = std::numeric_limits<float>::max(); // nano seconds
        GradientDescent gd;
    };
    
    struct PartitionningSpecs {
        PartitionningSpec & getWithSpread(bool spread) {
            if(spread) {
                return with_spread.avg_time_per_sample < without_spread.avg_time_per_sample ? with_spread : without_spread;
            }
            else {
                return without_spread;
            }
        }
        
        PartitionningSpec with_spread, without_spread;
    };
    
    template<typename T>
    struct PartitionAlgo {
        template<typename ...Args>
        static PartitionningSpecs run(Args... args) {
            return {};
        }
    };
    
    template<typename T>
    struct PartitionAlgo< PartitionnedFFTConvolution<T> > {
        static PartitionningSpecs run(int n_channels, int n_audio_cb_frames, int size_impulse_response) {
            assert(n_channels > 0);
            PartitionningSpecs res;
            {
                auto & spec = res.without_spread;
                spec.size = get_optimal_partition_size( spec.gd, n_channels, false, n_audio_cb_frames, size_impulse_response, spec.avg_time_per_sample );
            }
            
            if(n_channels > 1) {
                auto & spec = res.with_spread;
                spec.size = get_optimal_partition_size( spec.gd, n_channels, true, n_audio_cb_frames, size_impulse_response, spec.avg_time_per_sample );
            }
            
            return std::move(res);
        }
    };
    
}