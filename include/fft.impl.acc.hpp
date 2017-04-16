/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule {
    // implementation of Accelerate vDSP fft
    
    namespace accelerate {
        struct Tag {};
    }
    
    namespace fft {
        template<typename T>
        struct RealInput_<accelerate::Tag, T> {
            using type = std::vector<T>;
        };
        
        template<typename T>
        struct Context_<accelerate::Tag, T> {
            using type = accelerate::FFTSetup_<T>;
            
            static auto create(int lgSize) {
                return accelerate::API<T>::f_create_fftsetup(lgSize, kFFTRadix2);
            }
            static constexpr auto destroy = accelerate::API<T>::f_destroy_fftsetup;
        };

    }
    
    namespace accelerate {
        namespace fft {
            using namespace imajuscule::fft;

            template<typename T>
            using RealInput = typename RealInput_<Tag, T>::type;
            
            template<typename T>
            using Context = typename Context_<Tag, T>::type;
            
            template<typename T>
            using ScopedContext = ScopedContext_<Tag, T>;
        }
    }
}

