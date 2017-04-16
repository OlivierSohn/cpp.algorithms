/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule {
    namespace fft {
        
        template<typename T>
        static complex<T> make_root_of_unity(unsigned int index, unsigned int size) {
            using Tr = NumTraits<T>;
            return polar(-Tr::two() * static_cast<T>(M_PI) * index / size);
        }
        
        template<typename T>
        using FFTVec = typename std::vector<complex<T>>;

        template<typename T>
        using FFTIter = typename FFTVec<T>::iterator;
        
        template<typename T>
        void compute_roots_of_unity(unsigned int N, FFTVec<T> & res) {
            assert(is_power_of_two(N));
            auto n_roots = N/2;
            res.reserve(n_roots);
            for(unsigned int i=0; i<n_roots; ++i) {
                res.push_back(make_root_of_unity<T>(i,N));
            }
        }

        template<typename T>
        FFTVec<T> compute_roots_of_unity(unsigned int N) {
            FFTVec<T> res;
            compute_roots_of_unity(N, res);
            return std::move(res);
        }

        // https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm
        
        template<typename ROOTS_ITER, typename ITER1, typename ITER2>
        static inline void tukeyCooley(ROOTS_ITER root_it,
                         ITER1 it,
                         ITER2 result,
                         unsigned int N,
                         unsigned int stride) {
            if(N==0) {
                *result = *it;
                return;
            }
            auto double_stride = 2*stride;
            auto half_N = N/2;
            tukeyCooley(root_it, it         , result , half_N, double_stride );
            auto result2 = result + N;
            tukeyCooley(root_it, it + stride, result2, half_N, double_stride );
            
            for(unsigned int i=0;
                i<N;
                ++i, ++result, ++result2, root_it += stride)
            {
                auto & r1 = *result;
                auto & r2 = *result2;
                auto t = r1;
                r2 *= *root_it;
                r1 += r2;
                r2 = t - r2;
            }
        }

        template<typename T>
        struct Algo {
            using ROOTS_OF_UNITY = FFTVec<T>;
            using ROOTS_ITER = typename ROOTS_OF_UNITY::const_iterator;
            
            Algo() = default;
            
            Algo(ROOTS_OF_UNITY roots_of_unity) :
            roots_of_unity(std::move(roots_of_unity))
            {}
            
            void setRootsOfUnity(ROOTS_OF_UNITY roots) {
                roots_of_unity = std::move(roots);
            }
            
            template<typename ITER1, typename ITER2>
            void run(ITER1 it,
                     ITER2 result,
                     unsigned int N,
                     unsigned int stride) const {
                static_assert(std::is_same<T,typename ITER1::value_type::FPT>::value, "");
                static_assert(std::is_same<T,typename ITER2::value_type::FPT>::value, "");
                assert(N/2 == roots_of_unity.size());
                tukeyCooley(roots_of_unity.begin(),
                            it,
                            result,
                            N/2,
                            stride);
            }
            
        private:
            ROOTS_OF_UNITY roots_of_unity;
        };
        
        template<typename ITER>
        void compute_fft(unsigned int fft_length,
                         ITER signal_it,
                         ITER result_it) {
            using T = typename ITER::value_type;
            using FPT = typename T::FPT;
            auto roots = compute_roots_of_unity<FPT>(fft_length);
            Algo<FPT> a(std::move(roots));
            a.run(signal_it, result_it, fft_length, 1);
        }

        template<typename Iter>
        void normalize_fft(unsigned int fft_length,
                           Iter begin,
                           Iter end) {
            assert(fft_length);
            auto inv_l = 1. / fft_length;
            
            std::transform(begin, end,
                           begin, [inv_l](auto value) { return inv_l * value; });
        }
        
        template<typename T>
        struct FPT {
            using type = typename T::FPT;
        };
        
        template<>
        struct FPT<double> {
            using type = double;
        };

        template<>
        struct FPT<float> {
            using type = float;
        };
        
        template<typename ITER>
        void apply_hann_window(ITER it,
                               ITER end)
        {
            using V = typename ITER::value_type;
            using T = typename FPT<V>::type;
            auto NumTaps = std::distance(it, end);
            auto nyquist = NumTaps / 2;

            int i=0;
            for(; it != end; ++it) {
                // W(n) = cos(n/NumTaps · π/2)
                *it *= std::cos( std::abs(nyquist-i)/static_cast<T>(nyquist) * static_cast<T>(M_PI_2));
                ++i;
            }
        }

    } // NS fft
} // NS imajuscule

