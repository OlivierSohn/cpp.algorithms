/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule
{
    template<typename T>
    std::vector<complex<T>> complexify(std::vector<T> vec) {
        std::vector<complex<T>> ret;
        ret.reserve(vec.size());
        for(auto s : vec) {
            ret.push_back({s, 0});
        }
        return ret;
    }
    
    /*
     * cf. https://en.wikipedia.org/wiki/Overlap%E2%80%93add_method#math_Eq.1
     *
     * http://cnx.org/contents/tS7cr9DP@7/Fast-Convolution-Using-the-FFT
     *
     * todo use http://stackoverflow.com/questions/7972952/how-to-interpret-the-result-from-kissffts-kiss-fftr-fft-for-a-real-signal-fun
     * to optimize real ffts
     */
    template <typename T>
    struct FFTConvolution {
        FFTConvolution() {
            it = y.end();
        }
        
        bool empty() const { return fft_of_h.empty(); }
        
        void setCoefficients(std::vector<T> coeffs_) {
            N = coeffs_.size();
            auto N_nonzero_y = 2 * N - 1;
            
            auto fft_length = ceil_power_of_two(N_nonzero_y);
            
            fft_of_h.resize(fft_length, {0,0});
            fft_of_x.resize(fft_length, {0,0});
            result.resize(fft_length, {0,0});
            x.reserve(fft_length);
            {
                y.resize( 2 * N, {0,0});
                it = y.begin();
            }
            
            // pad impulse response with 0
            
            coeffs_.resize(fft_length, {});
            
            // compute fft of padded impulse response
            
            auto roots = fft::compute_roots_of_unity<double>(fft_length);
            fft.setRootsOfUnity(std::move(roots));
            auto coeffs = complexify(std::move(coeffs_));
            fft.run(coeffs.begin(), fft_of_h.begin(), fft_length, 1);
        }

        void step(T val) {
            x.emplace_back(val, 0);
            
            if(x.size() == N) {
                // pad x
                
                x.resize(fft_of_x.size(), {0,0});

                // do fft of x
                
                fft.run(x.begin(), fft_of_x.begin(), fft_of_x.size(), 1);
                
                // multiply fft_of_x by fft_of_h
                
                auto it_fft_x = fft_of_x.begin();
                auto end_fft_x = fft_of_x.end();

                auto it_fft_h = fft_of_h.begin();
                for(; it_fft_x != end_fft_x; ++it_fft_x, ++it_fft_h) {
                    auto & c = *it_fft_x;
                    c *= *it_fft_h;
                    // anticipating an operation required for inverse fft, cf. https://www.dsprelated.com/showarticle/800.php
                    c.convert_to_conjugate();
                }
                
                // https://www.dsprelated.com/showarticle/800.php
                fft.run(fft_of_x.begin(), result.begin(), fft_of_x.size(), 1);
                
                // in theory for inverse fft we should convert_to_conjugate the result
                // but it is supposed to be real numbers so the conjugation has no effect

                auto factor = 1 / static_cast<double>(fft_of_x.size());
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
            assert(it < y.begin() + N);
            assert(it >= y.begin());
            auto val = it->real();
            assert(std::abs(it->imag()) < 0.0001f);
            return val;
        }
        
    private:
        fft::Algo<double> fft;
        int N;
        
        std::vector<complex<double>> fft_of_h, x, fft_of_x, y, result; // todo use std::vector<T> when we have optimized real ffts
        decltype(y)::iterator it;
    };
}
