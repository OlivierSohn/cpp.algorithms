/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule {
    
    // implementation of custom (imajuscule) fft
    
    namespace imj {
        struct Tag {};
    }
    
    namespace fft {
        
        /*
         * Space complexity, for forward fft of real input of size N:
         *
         * input : 2*N
         */

        template<typename T>
        struct RealInput_<imj::Tag, T> {
            using type = std::vector<complex<T>>;
        };
    }
    
    namespace imj {
        namespace fft {
            using namespace imajuscule::fft;
            
            template<typename T>
            using RealInput = typename RealInput_<Tag, T>::type;
        }
    }
}

