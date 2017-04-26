/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule {
    
    // imajuscule's fft implementation
    
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
        struct RealSignal_<imj::Tag, T> {
            using type = a64::vector<complex<T>>;
            using iter = typename type::iterator;
            using const_iter = typename type::const_iterator;
            using value_type = typename type::value_type;

            static type make(a64::vector<T> reals) {
                type ret;
                ret.reserve(reals.size());
                for(auto r : reals) {
                    ret.emplace_back(r);
                }
                return std::move(ret);
            }
            
            static T get_signal(value_type const & c) {
                assert(std::abs(c.imag()) < 0.0001f);
                return c.real();
            }
            
            static void add_scalar_multiply(iter res, const_iter add1, const_iter add2, T const m, int N) {
                // res = m * (add1 + add2)
                
                for(int i=0; i<N; ++i, ++res, ++add1, ++add2) {
                    *res = m * (*add1 + *add2);
                }
            }
            
            static void copy(iter dest, const_iter from, int N) {
                memcpy(&*dest, &*from, N * sizeof(*dest));
            }

            static void zero(type & v) {
                std::fill(v.begin(), v.end(), value_type{});
            }
        };
        
        template<typename T>
        struct RealFBins_<imj::Tag, T> {
            using Tag = imj::Tag;
            using type = a64::vector<complex<T>>;
            
            static type make(type cplx) {
                return std::move(cplx);
            }

            static void mult_assign(type & v, type const & w) {
                auto it = v.begin();
                auto end = v.end();
                
                auto it_w = w.begin();
                for(; it != end; ++it, ++it_w) {
                    *it *= *it_w;
                }
            }
 
            static void zero(type & v) {
                complex<T> zero{};
                std::fill(v.begin(), v.end(), zero);
            }
            
            static void multiply_add(type & accum, type const & m1, type const & m2) {
                
                auto it_accum = accum.begin();
                
                auto it1 = m1.begin();
                auto end1 = m1.end();
                
                auto it2 = m2.begin();
                for(; it1 != end1; ++it2, ++it1, ++it_accum) {
                    assert(it_accum < accum.end());
                    *it_accum += *it1 * *it2;
                }
            }
            
            static std::pair<int, T> getMaxSquaredAmplitude(type const & v) {
                auto Max = T(0);
                
                int index = -1;
                int i=0;
                for( auto & c : v) {
                    auto M = norm(c);
                    if(M > Max) {
                        index = i;
                        Max = M;
                    }
                    ++i;
                }

                auto div = static_cast<T>(v.size()) * Algo_<Tag,T>::scale;

                return {index, Max/(div * div)};
            }
        };
        
        template<typename T>
        struct ImjContext {
            using vec_roots = a64::vector<complex<T>>;
            
            ImjContext() : roots(nullptr) {}
            ImjContext(vec_roots * roots) : roots(roots) {}
            
            operator bool() const {
                return !empty();
            }
            bool empty() const { return !roots; }
            void clear() { roots = nullptr; }
            
            vec_roots * getRoots() const { return roots; }
            vec_roots * editRoots() { return roots; }
        private:
            vec_roots * roots;
        };
        
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

        // https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm
        
        enum class FftType {
            FORWARD,
            INVERSE
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
            tukeyCooley<TYPE>(root_it, it         , result , half_N, double_stride );
            auto result2 = result + N;
            tukeyCooley<TYPE>(root_it, it + stride, result2, half_N, double_stride );
            
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
        struct Algo_<imj::Tag, T> {
            using RealInput  = typename RealSignal_ <imj::Tag, T>::type;
            using RealFBins  = typename RealFBins_<imj::Tag, T>::type;
            using Context    = typename Context_   <imj::Tag, T>::type;
            
            using Tr = NumTraits<T>;
            
            static constexpr auto scale = Tr::one();

            Algo_() = default;
            Algo_(Context c) : context(c) {}
            
            void setContext(Context c) {
                context = c;
            }

            void forward(RealInput const & input,
                         RealFBins & output,
                         unsigned int N) const
            {
                constexpr auto stride = 1;
                tukeyCooley<FftType::FORWARD>(context.getRoots()->begin(),
                                              input.begin(),
                                              output.begin(),
                                              N/2,
                                              stride);
            }
            
            void inverse(RealFBins const & input,
                         RealInput & output,
                         unsigned int N) const
            {
                constexpr auto stride = 1;
                tukeyCooley<FftType::INVERSE>(context.getRoots()->begin(),
                                              input.begin(),
                                              output.begin(),
                                              N/2,
                                              stride);
                // in theory for inverse fft we should convert_to_conjugate the result
                // but it is supposed to be real numbers so the conjugation would have no effect
                
#ifndef NDEBUG
                T M {};
                std::for_each(output.begin(), output.end(),
                              [&M](auto v) { M = std::max(M, std::abs(v.real())); } );
                for(auto const & r : output) {
                    if(M) {
                        assert(std::abs(r.imag()/M) < 1e-6);
                    }
                    else {
                        assert(std::abs(r.imag()) < 1e-6);
                    }
                }
#endif
            }
            
            Context context;
        };
        
        
        namespace slow_debug {
            
            template<typename CONTAINER>
            struct UnwrapFrequenciesRealFBins<imj::Tag, CONTAINER> {
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
            using RealInput = typename RealSignal_<Tag, T>::type;

            template<typename T>
            using RealFBins = typename RealFBins_<Tag, T>::type;
        
            template<typename T>
            using Context = typename Context_<Tag, T>::type;
            
            template<typename T>
            using ScopedContext = ScopedContext_<Tag, T>;
            
            template<typename T>
            using Algo = Algo_<Tag, T>;
        } // NS fft
    } // NS imj
} // NS imajuscule

