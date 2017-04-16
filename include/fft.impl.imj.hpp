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
        
        template<typename T>
        struct RealOutput_<imj::Tag, T> {
            using type = std::vector<complex<T>>;
        };
        
        template<typename T>
        struct ImjContext {
            using vec_roots = std::vector<complex<T>>;
            
            ImjContext(vec_roots * roots) : roots(roots) {}
            
            vec_roots * getRoots() const { return roots; }
            vec_roots * editRoots() { return roots; }
        private:
            vec_roots * roots;
        };
        
        template<typename T>
        static complex<T> make_root_of_unity(unsigned int index, unsigned int size) {
            using Tr = NumTraits<T>;
            return polar(-Tr::two() * static_cast<T>(M_PI) * index / size);
        }
        
        template<typename T>
        void compute_roots_of_unity(unsigned int N, std::vector<complex<T>> & res) {
            assert(is_power_of_two(N));
            auto n_roots = N/2;
            res.reserve(n_roots);
            for(unsigned int i=0; i<n_roots; ++i) {
                res.push_back(make_root_of_unity<T>(i,N));
            }
        }
        
        template<typename T>
        auto compute_roots_of_unity(unsigned int N) {
            std::vector<complex<T>> res;
            compute_roots_of_unity(N, res);
            return std::move(res);
        }
        
        // https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm
        
        enum class FftType {
            FORWARD,
            INVERSE,
            PRIVATE_INVERSE_SUBTREE
        };
        template<FftType TYPE, typename ROOTS_ITER, typename ITER1, typename ITER2>
        static inline void tukeyCooley(ROOTS_ITER root_it,
                                       ITER1 it,
                                       ITER2 result,
                                       unsigned int N,
                                       unsigned int stride) {
            if(N==0) {
                if(TYPE == FftType::FORWARD) {
                    *result = *it;
                }
                else {
                    // conjugate
                    auto & c = *it;
                    *result = { c.real(), -c.imag() };
                }
                return;
            }
            auto double_stride = 2*stride;
            auto half_N = N/2;
            if(TYPE == FftType::PRIVATE_INVERSE_SUBTREE) {
                tukeyCooley<FftType::PRIVATE_INVERSE_SUBTREE>
                (root_it, it         , result , half_N, double_stride );
            }
            else {
                tukeyCooley<FftType::FORWARD>
                (root_it, it         , result , half_N, double_stride );
            }
            auto result2 = result + N;
            if(TYPE == FftType::FORWARD) {
                tukeyCooley<FftType::FORWARD>
                (root_it, it + stride, result2, half_N, double_stride );
            }
            else {
                tukeyCooley<FftType::PRIVATE_INVERSE_SUBTREE>
                (root_it, it + stride, result2, half_N, double_stride );
            }
            
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
        struct Context_<imj::Tag, T> {
            using type = ImjContext<T>;
            using InnerCtxt = typename type::vec_roots;
            
            static auto create(int size) {
                auto pv = new InnerCtxt();
                compute_roots_of_unity(size, *pv);
                return type(pv);
            }

            static void destroy(type c) {
                delete c.editRoots();
            }
        };
        
        template<typename T>
        struct Algo_<imj::Tag, T> {
            using RealInput  = typename RealInput_ <imj::Tag, T>::type;
            using RealOutput = typename RealOutput_<imj::Tag, T>::type;
            using Context    = typename Context_   <imj::Tag, T>::type;
            
            using Tr = NumTraits<T>;
            
            static constexpr auto scale = Tr::one();

            Algo_(Context c) : ctxt(c) {}
            
            void forward(RealInput const & input,
                         RealOutput & output,
                         unsigned int N) const
            {
                constexpr auto stride = 1;
                tukeyCooley<FftType::FORWARD>(ctxt.getRoots()->begin(),
                                              input.begin(),
                                              output.begin(),
                                              N/2,
                                              stride);
            }
            void inverse(RealInput const & input,
                         RealOutput & output,
                         unsigned int N) const
            {
                constexpr auto stride = 1;
                tukeyCooley<FftType::INVERSE>(ctxt.getRoots()->begin(),
                                              input.begin(),
                                              output.begin(),
                                              N/2,
                                              stride);
            }
            
            Context ctxt;
        };
        
        
        namespace slow_debug {
            
            template<typename CONTAINER>
            struct UnwrapFrequencies<imj::Tag, CONTAINER> {
                static auto run(CONTAINER container, int N) {
                    return std::move(container);
                }
            };
            
            template<typename CONTAINER>
            struct UnwrapSignal<imj::Tag, CONTAINER> {
                static auto run(CONTAINER container, int N) {
                    return std::move(container);
                }
            };
            
        } // NS slow_debug
    } // NS fft
    
    namespace imj {
        namespace fft {
            using namespace imajuscule::fft;
            
            // this part could be #included to avoid repetitions
            
            template<typename T>
            using RealInput = typename RealInput_<Tag, T>::type;

            template<typename T>
            using RealOutput = typename RealOutput_<Tag, T>::type;
        
            template<typename T>
            using Context = typename Context_<Tag, T>::type;
            
            template<typename T>
            using ScopedContext = ScopedContext_<Tag, T>;
            
            template<typename T>
            using Algo = Algo_<Tag, T>;
        } // NS fft
    } // NS imj
} // NS imajuscule

