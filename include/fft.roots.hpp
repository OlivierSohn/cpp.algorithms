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
        void compute_roots_of_unity(unsigned int N, a64::vector<complex<T>> & res) {
            assert(is_power_of_two(N));
            auto n_roots = N/2;
            res.reserve(n_roots);
            for(unsigned int i=0; i<n_roots; ++i) {
                res.push_back(make_root_of_unity<T>(i,N));
            }
        }
        
        template<typename T>
        auto compute_roots_of_unity(unsigned int N) {
            a64::vector<complex<T>> res;
            compute_roots_of_unity(N, res);
            return std::move(res);
        }
        
    } // NS fft
    
} // NS imajuscule

