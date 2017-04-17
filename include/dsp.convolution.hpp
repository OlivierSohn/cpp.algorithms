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
    template <typename Parent, typename Tag>
    struct FFTConvolutionBase : public Parent {
        using T = typename Parent::FPT;
        using FPT = T;
        using FFT = typename Parent::FFT;
        using FFTAlgo = typename fft::Algo_<Tag, FPT>;
        using Contexts = fft::Contexts_<Tag, FPT>;
        
        using Parent::get_fft_length;
        using Parent::getComputationPeriodicity;
        using Parent::doSetCoefficients;
        using Parent::compute_convolution;
        
        FFTConvolutionBase() {
            it = y.end();
        }
        
        void setCoefficients(std::vector<T> coeffs_) {
            
            auto const N = coeffs_.size();
            auto const fft_length = get_fft_length(N);
            fft.setContext(Contexts::getInstance().getBySize(fft_length));
            
            result.resize(fft_length, {0,0});
            x.reserve(fft_length);
            {
                y.resize( fft_length, {0,0});
                it = y.begin();
            }
            
            doSetCoefficients(fft, std::move(coeffs_));
        }
        
        void step(T val) {
            x.emplace_back(val, 0);
            
            auto N = getComputationPeriodicity();
            if(unlikely(x.size() == N)) {
                // pad x
                
                x.resize(get_fft_length(), {0,0});
                
                auto const & frequencies = compute_convolution(fft, x);
                
                fft.inverse(frequencies, result, get_fft_length());
                
                // in theory for inverse fft we should convert_to_conjugate the result
                // but it is supposed to be real numbers so the conjugation would have no effect
                
#ifndef NDEBUG
                for(auto const & r : result) {
                    assert(r.imag() < 0.01f);
                }
#endif

                auto factor = 1 / (FFTAlgo::scale * static_cast<T>(get_fft_length()));
                // for efficiency reasons we will scale the result in the following loop to avoid an additional traversal
                
                // y = mix first part of result with second part of previous result
                
                auto it_res = result.begin();
                auto it_y = y.begin();
                auto end_y = it_y + N;
                auto it_y_prev = end_y;
                auto end_y_prev = it_y_prev + N;
                
                for(; it_y != end_y; ++it_y, ++it_y_prev, ++it_res) {
                    *it_y = factor * (*it_y_prev + *it_res);
                }
                
                // store second part of result for later
                for(; it_y != end_y_prev; ++it_y, ++it_res) {
                    *it_y = *it_res;
                }
                
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
            auto val = it->real();
            assert(std::abs(it->imag()) < 0.0001f);
            return val;
        }
        
    private:
        FFTAlgo fft;
        FFT x, y, result;
        typename decltype(y)::iterator it;
    };
    
    template <typename T>
    struct FFTConvolutionCRTP {
        using FPT = T;
        using FFT = typename fft::FFTVec<T>;

        auto getComputationPeriodicity() const { return N; }
        
        bool empty() const { return fft_of_h.empty(); }
        
        auto get_fft_length(int N) {
            auto N_nonzero_y = 2 * N - 1;
            return ceil_power_of_two(N_nonzero_y);
        }
        
        auto get_fft_length() {
            return fft_of_x.size();
        }
        
        template<typename FFTAlgo>
        void doSetCoefficients(FFTAlgo const & fft, std::vector<T> coeffs_) {
            N = coeffs_.size();
            
            auto fft_length = get_fft_length(N);
            
            fft_of_h.resize(fft_length, {0,0});
            fft_of_x.resize(fft_length, {0,0});
            
            // pad impulse response with 0
            
            coeffs_.resize(fft_length, {});
            
            // compute fft of padded impulse response
            
            auto coeffs = complexify<FPT>(coeffs_.begin(), coeffs_.end());
            fft.forward(coeffs, fft_of_h, fft_length);
        }
        
    protected:

        template<typename FFTAlgo>
        auto const & compute_convolution(FFTAlgo const & fft, FFT const & x) {
            // do fft of x
            
            fft.forward(x, fft_of_x, get_fft_length());
            
            // multiply fft_of_x by fft_of_h
            
            auto it_fft_x = fft_of_x.begin();
            auto end_fft_x = fft_of_x.end();
            
            auto it_fft_h = fft_of_h.begin();
            for(; it_fft_x != end_fft_x; ++it_fft_x, ++it_fft_h) {
                *it_fft_x *= *it_fft_h;
            }
            
            return fft_of_x;
        }
      
    private:
        int N;
        
        FFT fft_of_h, fft_of_x; // todo use std::vector<T> when we have optimized real ffts
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

    template <typename T>
    struct PartitionnedFFTConvolutionCRTP {
        using FPT = T;
        using FFT = typename fft::FFTVec<T>;
        
        auto get_fft_length() { assert(partition_size > 0); return 2 * partition_size; }
        auto get_fft_length(int) { return get_fft_length(); }
        
        bool empty() const { return ffts_of_partitionned_h.empty(); }

        auto getComputationPeriodicity() { return partition_size; }
        
        void set_partition_size(int sz) {
            assert(sz > 0);
            partition_size = sz;
            assert(is_power_of_two(sz));
        }

        template<typename FFTAlgo>
        void doSetCoefficients(FFTAlgo const & fft, std::vector<T> coeffs_) {
            
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
                fft_of_partitionned_h.resize(fft_length, {0,0});
            }
            for(auto & fft_of_delayed_x : ffts_of_delayed_x) {
                fft_of_delayed_x.resize(fft_length, {0,0});
            }
            
            work.resize(fft_length);
            
            // compute fft of padded impulse response
            
            auto it_coeffs = coeffs_.begin();
            for(auto & fft_of_partitionned_h : ffts_of_partitionned_h) {
                auto end_coeffs = it_coeffs + partition_size;
                assert(end_coeffs <= coeffs_.end());
                auto coeffs = complexify<FPT>(it_coeffs, end_coeffs);
                it_coeffs = end_coeffs;

                // pad partitionned impulse response with 0
                
                coeffs.resize(fft_length, {0,0});
                fft.forward(coeffs, fft_of_partitionned_h, fft_length);
            }
            assert(it_coeffs == coeffs_.end());
        }

    protected:

        template<typename FFTAlgo>
        auto const & compute_convolution(FFTAlgo const & fft, FFT const & x)
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
            
            std::fill(work.begin(), work.end(), complex<T>(0,0));
            
            ffts_of_delayed_x.for_each_bkwd( [this, &it_fft_of_partitionned_h] (auto const & fft_of_delayed_x) {
                assert(it_fft_of_partitionned_h < ffts_of_partitionned_h.end());
                
                auto const & fft_of_partitionned_h = *it_fft_of_partitionned_h;
                
                // multiply fft_of_delayed_x by fft_of_partitionned_h
                
                auto it_work = work.begin();
                
                auto it_fft_x = fft_of_delayed_x.begin();
                auto end_fft_x = fft_of_delayed_x.end();
                
                auto it_fft_h = fft_of_partitionned_h.begin();
                for(; it_fft_x != end_fft_x; ++it_fft_x, ++it_fft_h, ++it_work) {
                    assert(it_work < work.end());
                    *it_work += *it_fft_x * *it_fft_h;
                }
                
                ++ it_fft_of_partitionned_h;
            });
            
            return work;
        }

    private:
        int partition_size = -1;
        cyclic<FFT> ffts_of_delayed_x;
        std::vector<FFT> ffts_of_partitionned_h;
        
        FFT work;
    };

    using FFTAlgoTag = imj::Tag;
    
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
    template <typename T>
    using FFTConvolution = FFTConvolutionBase< FFTConvolutionCRTP<T> , FFTAlgoTag >;

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
    template <typename T>
    using PartitionnedFFTConvolution = FFTConvolutionBase< PartitionnedFFTConvolutionCRTP<T>, FFTAlgoTag >;
    
    
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
