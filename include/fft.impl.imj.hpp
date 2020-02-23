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
            
            /* we need complex<T> (type) to do the ftt and ifft, but after that,
             we can use T (outputType) */
            
            using outputType = a64::vector<T>;
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
            
            template<typename T2>
            static void add_assign_output(T2 * __restrict res,
                                          T const * __restrict add,
                                          int const N) {
                for(auto const * resEnd = res + N;
                    res != resEnd;
                    ++res, ++add)
                {
                    *res += *add;
                }
            }
            
            static void add_assign(T * __restrict res,
                                   value_type const * __restrict add,
                                   int const start,
                                   int const N) {
                add += start;
                for(auto const * resEnd = res + N;
                    res != resEnd;
                    ++res, ++add)
                {
                    *res += add->real();
                }
            }

            static void copy(value_type * __restrict dest,
                             value_type const * __restrict from,
                             int N) {
                // TODO optimize ?
                memcpy(dest,
                       from,
                       N * sizeof(value_type));
            }

            template<typename T2>
            static void copyOutputToOutput(T2 * __restrict dest,
                                           T const * __restrict from,
                                           int N) {
                // TODO optimize ?
                if constexpr(std::is_same_v<T, T2>) {
                    memcpy(dest,
                           from,
                           N * sizeof(T));
                }
                else {
                    for(int i=0; i!=N; ++i) {
                        dest[i] = from[i];
                    }
                }
            }

            template<typename TSource>
            static void copyFromInput(value_type * __restrict dest,
                                      TSource const * __restrict from,
                                      int N) {
                for(int i=0; i!=N; ++i) {
                    dest[i] = value_type(from[i]);
                }
            }
            
            static void copyToOutput(T * __restrict dest,
                                     value_type const * __restrict from,
                                     int const start,
                                     int const N) {
                from += start;
                for(int i=0; i!=N; ++i) {
                    dest[i] = from[i].real();
                }
            }

            static void zero_n_raw_output(T * p, int n) {
                std::fill(p, p+n, T{});
            }
            static void zero_n_raw(complex<T> * p, int n) {
                std::fill(p, p+n, value_type{});
            }
            static void zero_n(type & v, int n) {
                zero_n_raw(v.data(), n);
            }
            static void zero(type & v) {
                zero_n(v, v.size());
            }
            
            static void dotpr(complex<T> const * const a,
                              T const * const b,
                              T * res,
                              int n) {
                T r{};
                for(int i=0; i<n; ++i) {
                    Assert(0==a->imag());
                    r += a[i].real() * b[i];
                }
                assert(res);
                *res = r;
            }
        };

        template<typename T, template<typename> typename Allocator>
        struct RealFBins_<imj::Tag, T, Allocator> {
            using Tag = imj::Tag;
            using type = std::vector<complex<T>, Allocator<complex<T>>>;

            static type make(type cplx) {
                return std::move(cplx);
            }

          static void scale(type & v, T const s) {
            auto * __restrict it = v.begin().base();
            auto * __restrict end = v.end().base();
            
            for(; it != end; ++it) {
              *it *= s;
            }
          }

            static void mult_assign(type & v, type const & w) {
                auto * __restrict it = v.begin().base();
                auto * __restrict end = v.end().base();

                auto * __restrict it_w = w.begin().base();
                for(; it != end; ++it, ++it_w) {
                    *it *= *it_w;
                }
            }

            static void zero(type & v) {
                complex<T> zero{};
                std::fill(v.begin(), v.end(), zero);
            }
            
            static void multiply(complex<T> * const __restrict res,
                                 complex<T> const * const __restrict m1,
                                 complex<T> const * const __restrict m2,
                                 int N) {
                Assert(N > 0);
                for(int i=0, end = 2*N; i != end; ++i)
                {
                    res[i] = m1[i] * m2[i];
                }
            }
            static void multiply_add(complex<T> * const __restrict res,
                                     complex<T> const * const __restrict m1,
                                     complex<T> const * const __restrict m2,
                                     int N) {
                Assert(N > 0);
                for(int i=0, end = 2*N; i != end; ++i)
                {
                    res[i] += m1[i] * m2[i];
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
            
            static void setFirstReal(type & v, T value) {
                if(v.empty()) {
                    throw std::logic_error("setFirstReal on empty");
                }
                v[0] = complex{value};
            }
            static T getFirstReal(type const & v) {
                if(v.empty()) {
                    throw std::logic_error("getFirstReal on empty");
                }
                return v[0].real();
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

      template<FftType TYPE, typename T>
      struct TukeyCooley {
        complex<T> * const root;
        
        // N is 'result' size
        void run(complex<T> const * const __restrict input,
                 complex<T> * __restrict result,
                 unsigned int const N) const {
          tukeyCooley(input, result, N/2, 1);
        }
      private:
    
      /*       
       Memory acces pattern:
       
       r1/i1 r2/ir2 .... r2n/i2n
       
       accessed in this order:
       
       when N=0
       
       r1/i1 rn/in
       r(n/2)/i(n/2) r(3*n/2)/i(3*n/2)
       
       r(n/4)
       ...
       */
        void tukeyCooley(complex<T> const * const __restrict it,
                         complex<T> * __restrict result,
                         unsigned int const N,
                         unsigned int const stride) const {
          if(N==0) {
            if constexpr (TYPE == FftType::FORWARD) {
              *result = *it;
            }
            else {
              *result = conj(*it);
            }
            return;
          }
          auto const double_stride = 2*stride;
          auto const half_N = N/2;
          // computes first half of result
          // using input with offset 0
          tukeyCooley(it         , result , half_N, double_stride );
          auto * __restrict result2 = result + N;
          // computes second half of result
          // using input with offset stride
          tukeyCooley(it + stride, result2, half_N, double_stride );
          
          // full result by mixing the 2 halves
          complex<T> * __restrict root_it = root;
          for(;result != result2;
              ++result, root_it += stride)
          {
            auto const t = result[N] * *root_it;
            result[N] = result[0] - t;
            result[0] += t;
          }
        }
      };

        template<typename T>
        struct Algo_<imj::Tag, T> {
            static constexpr bool inplace_dft = false;

            using FPT = T;

            using RealInput  = typename RealSignal_ <imj::Tag, T>::type;
            using RealValue = typename RealInput::value_type;
            
            template<template<typename> typename Allocator>
            using RealFBins  = typename RealFBins_<imj::Tag, T, Allocator>::type;

            using Context    = typename Context_   <imj::Tag, T>::type;

            using Tr = NumTraits<T>;

            static constexpr auto scale = Tr::one();

            Algo_() = default;
            Algo_(Context c) : context(c) {}

            void setContext(Context c) {
                context = c;
            }

          void forward(typename RealInput::const_iterator inputBegin,
                         complex<T> * output,
                         unsigned int N) const
          {
            auto * const rootPtr = context.getRoots()->begin().base();
            TukeyCooley<FftType::FORWARD, T>
            algo{rootPtr};
            
            algo.run(inputBegin.base(),
                     output,
                     N);
          }

            void inverse(complex<T> const * input,
                         RealValue * output,
                         unsigned int N) const
            {
              auto * const rootPtr = context.getRoots()->begin().base();
              TukeyCooley<FftType::INVERSE, T>
              algo{rootPtr};
              
              algo.run(input,
                       output,
                       N);

                // in theory for inverse fft we should convert_to_conjugate the result
                // but it is supposed to be real numbers so the conjugation would have no effect

#ifndef NDEBUG
                T M {};
                for(int i=0; i<N; ++i) {
                    M = std::max(M, std::abs(output[i].real()));
                }
                constexpr auto epsilon = 1000 * std::numeric_limits<T>::epsilon();
                for(int i=0; i<N; ++i) {
                    if(M) {
                        assert(std::abs(output[i].imag()/M) < epsilon);
                    }
                    else {
                        assert(std::abs(output[i].imag()) < epsilon);
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

            template<typename T, template<typename> typename Allocator>
            using RealFBins = typename RealFBins_<Tag, T, Allocator>::type;

            template<typename T>
            using Context = typename Context_<Tag, T>::type;

            template<typename T>
            using ScopedContext = ScopedContext_<Tag, T>;

            template<typename T>
            using Algo = Algo_<Tag, T>;
        } // NS fft
    } // NS imj
} // NS imajuscule
