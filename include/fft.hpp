/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule {
    namespace fft {
        
        template<typename T>
        using FFTVec = typename std::vector<complex<T>>;

        template<typename T>
        using FFTIter = typename FFTVec<T>::iterator;
        
        // TODO replace by others :
        // this does forward fft only so we wil need to remove the conjugations and use inverse where needed.
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
                tukeyCooley<FftType::FORWARD>(roots_of_unity.begin(),
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
            
            using namespace imj::fft;
            
            ScopedContext<FPT> scoped_context(fft_length);
            Algo<FPT> fft(scoped_context.get());
            fft.run(signal_it, result_it, fft_length);
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

