/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule
{
    template<typename T, typename ITER = typename std::vector<T>::iterator>
    std::vector<complex<T>> complexify(ITER it, ITER end) {
        std::vector<complex<T>> ret;
        ret.reserve(std::distance(it, end));
        for(; it!=end; ++it) {
            ret.push_back({*it, 0});
        }
        return ret;
    }
    
    /*
     * cf. https://en.wikipedia.org/wiki/Overlap%E2%80%93add_method#math_Eq.1
     */
    template <typename Parent>
    struct FFTConvolutionBase : public Parent {
        using T = typename Parent::FPT;
        using FPT = T;
        using FFT = typename Parent::FFT;
        using FFTAlgo = typename fft::Algo<FPT>;
        
        using Parent::get_fft_length;
        using Parent::getComputationPeriodicity;
        using Parent::doSetCoefficients;
        using Parent::compute_convolution;
        
        FFTConvolutionBase() {
            it = y.end();
        }
        
        void setCoefficients(std::vector<T> coeffs_) {
            
            auto N = coeffs_.size();
            auto fft_length = get_fft_length(N);
            auto roots = fft::compute_roots_of_unity<T>(fft_length);
            fft.setRootsOfUnity(std::move(roots));
            
            result.resize(fft_length, {0,0});
            x.reserve(fft_length);
            {
                y.resize( 2 * N, {0,0});
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
                
                auto it_ = compute_convolution(fft, x);
                
                // finish inverse fft
                
                fft.run(it_, result.begin(), get_fft_length(), 1);
                
                // in theory for inverse fft we should convert_to_conjugate the result
                // but it is supposed to be real numbers so the conjugation would have no effect
                
#ifndef NDEBUG
                for(auto const & r : result) {
                    assert(r.imag() < 0.000001f);
                }
#endif

                auto factor = 1 / static_cast<T>(get_fft_length());
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
        
        void doSetCoefficients(fft::Algo<FPT> const & fft, std::vector<T> coeffs_) {
            N = coeffs_.size();
            
            auto fft_length = get_fft_length(N);
            
            fft_of_h.resize(fft_length, {0,0});
            fft_of_x.resize(fft_length, {0,0});
            
            // pad impulse response with 0
            
            coeffs_.resize(fft_length, {});
            
            // compute fft of padded impulse response
            
            auto coeffs = complexify<FPT>(coeffs_.begin(), coeffs_.end());
            fft.run(coeffs.begin(), fft_of_h.begin(), fft_length, 1);
        }
        
        
    protected:
        auto compute_convolution(fft::Algo<FPT> const & fft, FFT const & x) {
            // do fft of x
            
            fft.run(x.begin(), fft_of_x.begin(), get_fft_length(), 1);
            
            // multiply fft_of_x by fft_of_h
            
            auto it_fft_x = fft_of_x.begin();
            auto end_fft_x = fft_of_x.end();
            
            auto it_fft_h = fft_of_h.begin();
            for(; it_fft_x != end_fft_x; ++it_fft_x, ++it_fft_h) {
                auto & c = *it_fft_x;
                c *= *it_fft_h;
                
                // start inverse fft (anticipation step here, cf. https://www.dsprelated.com/showarticle/800.php)
                c.convert_to_conjugate();
            }
            
            return fft_of_x.begin();
        }
      
    private:
        int N;
        
        FFT fft_of_h, fft_of_x; // todo use std::vector<T> when we have optimized real ffts
    };
    
    // we compute a table of "best lg2_partition_size" depending on:
    // - number of channels (1 and 2 for now)
    // - the number of frames computed per audio callback (32 and up)
    // - the ceil_power_2(length of impulse)
    //
    // we use gradient descent to find the best size.
    //
    // we redo the computation in mode "N_Channel * A < PART_L", knowing that in that case,
    // we have an advantage because we can distribute the various channels across successive callback calls
    // provided the "phases" of the different algorithms are well-spaced.
    // So in that mode, we divide the time by n_channels.
    //
    // We just need to measure one computation with ffts and time it in worst case, when nothing is in the cache.
    // To do that we need to explicitely pollute the cache.
    //
    // when deducing the right partition size, we don't interpolate on number of frames,
    // we take the best result between :
    //    ceil power of 2 in mode "N_Channel * A < PART_L"
    //    floor power of 2 in normal mode

    /*
     * Partitionned convolution, cf. http://www.ericbattenberg.com/school/partconvDAFx2011.pdf
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
    template <typename T, int LG2_PARTITION_SIZE>
    struct PartitionnedFFTConvolutionCRTP {
        using FPT = T;
        using FFT = typename fft::FFTVec<T>;
        
        static_assert(LG2_PARTITION_SIZE >= 0, "");
        
        static constexpr auto PARTITION_SIZE = pow2(LG2_PARTITION_SIZE);
        
        static constexpr auto fft_length = 2 * PARTITION_SIZE;
        
        constexpr auto get_fft_length() { return fft_length; }
        constexpr auto get_fft_length(int) { return get_fft_length(); }
        
        bool empty() const { return ffts_of_partitionned_h.empty(); }

        static auto getComputationPeriodicity() { return PARTITION_SIZE; }
        
        void doSetCoefficients(fft::Algo<FPT> const & fft, std::vector<T> coeffs_) {
            
            auto N = coeffs_.size();
            auto n_partitions = N/PARTITION_SIZE;
            if(n_partitions * PARTITION_SIZE != N) {
                // one partition is partial...
                assert(n_partitions * PARTITION_SIZE < N);
                ++n_partitions;
                // ... pad it with zeros
                coeffs_.resize(n_partitions * PARTITION_SIZE, {});
            }
            ffts_of_delayed_x.resize(n_partitions);
            ffts_of_partitionned_h.resize(n_partitions);
            
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
                auto end_coeffs = it_coeffs + PARTITION_SIZE;
                assert(end_coeffs <= coeffs_.end());
                auto coeffs = complexify<FPT>(it_coeffs, end_coeffs);
                it_coeffs = end_coeffs;

                // pad partitionned impulse response with 0
                
                coeffs.resize(fft_length, {0,0});
                fft.run(coeffs.begin(), fft_of_partitionned_h.begin(), fft_length, 1);
            }
            assert(it_coeffs == coeffs_.end());
        }

    protected:
        auto compute_convolution(fft::Algo<FPT> const & fft, FFT const & x) {
            
            // do fft of x
            {
                auto & oldest_fft_of_delayed_x = *ffts_of_delayed_x.cycleEnd();
                assert(fft_length == oldest_fft_of_delayed_x.size());
                fft.run(x.begin(), oldest_fft_of_delayed_x.begin(), fft_length, 1);
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
            
            // start inverse fft (https://www.dsprelated.com/showarticle/800.php)
            
            for(auto & c : work) {
                c.convert_to_conjugate();
            }
            
            return work.begin();
        }

    private:
        cyclic<FFT> ffts_of_delayed_x;
        std::vector<FFT> ffts_of_partitionned_h;
        
        FFT work;
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
    template <typename T>
    using FFTConvolution = FFTConvolutionBase< FFTConvolutionCRTP<T> >;

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
     *                also we could optimize more globally, taking into account that we have one reverb per channel:
     *                when N_Channel * A < PART_L we can distribute the computes over the different callbacks calls, provided the "phases"
     *                of the different algorithms are well-spaced.
     *
     */
    template <typename T, int LG2_PARTITION_SIZE>
    using PartitionnedFFTConvolution = FFTConvolutionBase< PartitionnedFFTConvolutionCRTP<T, LG2_PARTITION_SIZE> >;
    
}
