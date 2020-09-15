/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule {

namespace imj2 {
struct Tag {};
}

static inline unsigned int bitReverse(unsigned int b,
                               unsigned int log2N)
{
    b = (((b & 0xaaaaaaaa) >> 1) | ((b & 0x55555555) << 1));
    b = (((b & 0xcccccccc) >> 2) | ((b & 0x33333333) << 2));
    b = (((b & 0xf0f0f0f0) >> 4) | ((b & 0x0f0f0f0f) << 4));
    b = (((b & 0xff00ff00) >> 8) | ((b & 0x00ff00ff) << 8));
    b = ((b >> 16) | (b << 16)) >> (32 - log2N);
    return b;
}

    namespace fft {

        template<typename T>
        struct RealSignal_<imj2::Tag, T> {

            using type = a64::vector<T>;
            using outputType = type;
            using iter = typename type::iterator;
            using const_iter = typename type::const_iterator;
            using value_type = typename type::value_type;
            
            static type make(a64::vector<T> reals) {
                return reals;
            }

            static T get_signal(value_type const & c) {
                return c;
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
                                   complex<T> const * __restrict add,
                                   unsigned int const addSz,
                                   unsigned int const start,
                                   unsigned int const N) {
                unsigned int const log2addSz = [addSz]() {
                    assert(addSz);
                    if constexpr (overlapMode == Overlap::Add) {
                        return power_of_two_exponent(addSz);
                    }
                    else {
                        return 1 + power_of_two_exponent(addSz);
                    }
                }();
                
                for(unsigned int i=start, end=start+N; i!=end; ++i, ++res) {
                    unsigned int b = bitReverse(i, log2addSz);

                    assert(std::abs(add[b].imag()) < 0.0001f);

                    *res += add[b].real();
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
                    dest[i] = from[i];
                }
            }
            
            static void copyToOutput(T * __restrict dest,
                                     complex<T> const * __restrict from,
                                     unsigned int const fromSz,
                                     unsigned int const start,
                                     unsigned int const N) {
                unsigned int const log2fromSz = [fromSz]() {
                    assert(fromSz);
                    if constexpr (overlapMode == Overlap::Add) {
                        return power_of_two_exponent(fromSz);
                    }
                    else {
                        return 1 + power_of_two_exponent(fromSz);
                    }
                }();
                
                for(unsigned int i=start, end=start+N; i!=end; ++i, ++dest) {
                    unsigned int b = bitReverse(i, log2fromSz);

                    assert(std::abs(from[b].imag()) < 0.0001f);

                    *dest = from[b].real();
                }
            }

            static void zero_n_raw_output(T * p, int n) {
                std::fill(p, p+n, T{});
            }
            static void zero_n_raw(T * p, int n) {
                std::fill(p, p+n, T{});
            }
            static void zero_n(type & v, int n) {
                zero_n_raw(v.data(), n);
            }
            static void zero(type & v) {
                zero_n(v, v.size());
            }
            
            static void dotpr(T const * const a,
                              T const * const b,
                              T * res,
                              int n) {
                T r{};
                for(int i=0; i<n; ++i) {
                    r += a[i] * b[i];
                }
                assert(res);
                *res = r;
            }
        };

        template<typename T, template<typename> typename Allocator>
        struct RealFBins_<imj2::Tag, T, Allocator> {
            using Tag = imj2::Tag;
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
        struct Context_<imj2::Tag, T> {
            using type = int;

            static auto create(int size) {
                return int(0);
            }

            static void destroy(type c) {
            }
        };

        template<typename T>
        struct Algo_<imj2::Tag, T> {
            static constexpr bool inplace_dft = true;
            
            using FPT = T;

            using RealInput  = typename RealSignal_ <imj2::Tag, T>::type;
            using RealValue = typename RealInput::value_type;
            
            template<template<typename> typename Allocator>
            using RealFBins  = typename RealFBins_<imj2::Tag, T, Allocator>::type;

            using Context    = typename Context_   <imj2::Tag, T>::type;

            using Tr = NumTraits<T>;

            static constexpr auto scale = Tr::one();

            Algo_() = default;
            Algo_(Context c) {}

            void setContext(Context c) const {
            }
            
            void forward(typename RealInput::const_iterator inputBegin,
                         complex<T> * x,
                         unsigned int const N) const
            {
                // copy input to output
                
                for(int i=0; i<N; ++inputBegin, ++i) {
                    x[i] = complex<T>{*inputBegin};
                }
                
                dft_inplace(x, N);
                bit_reverse(x, N);
            }
            void inverse(complex<T> * x,
                         unsigned int N) const
            {
                // conjugate input
                
                for(int i=0; i<N; ++i) {
                    x[i].convert_to_conjugate();
                }

                dft_inplace(x, N);
            }
            
            static T extractRealOutput(complex<T> * x,
                                       unsigned int b,
                                       unsigned int N)
            {
                unsigned int log2N = power_of_two_exponent(N);

                b = bitReverse(b, log2N);

                return x[b].real();
            }
            
        private:
            void dft_inplace(complex<T> * x,
                             unsigned int const N) const
            {
                // https://rosettacode.org/wiki/Fast_Fourier_transform#C.2B.2B
                
                // DFT
                unsigned int k = N, n;
                double thetaT = 3.14159265358979323846264338328L / N;
                complex<T> phiT = complex<T>(cos(thetaT), -sin(thetaT)), T2;
                while (k > 1)
                {
                    n = k;
                    k >>= 1;
                    phiT = phiT * phiT;
                    T2 = complex<T>{1.0};
                    for (unsigned int l = 0; l < k; l++)
                    {
                        for (unsigned int a = l; a < N; a += n)
                        {
                            unsigned int b = a + k;
                            complex<T> t = x[a] - x[b];
                            x[a] += x[b];
                            x[b] = t * T2;
                        }
                        T2 *= phiT;
                    }
                }
            }
            
            template<typename TT>
            void bit_reverse(TT * x, unsigned int const N) const
            {
                // Reverse bits
                unsigned int m = power_of_two_exponent(N);
                for(unsigned int i=0; i<N; ++i) {
                    unsigned int b = bitReverse(i, m);
                    
                    if(b<i) {
                        std::swap(x[i], x[b]);
                    }
                }
            }
        };

        template<typename CONTAINER>
        struct UnwrapFrequenciesRealFBins<imj2::Tag, CONTAINER> {
            static auto run(CONTAINER container, int N) {
                return container;
            }
        };

        template<typename CONTAINER>
        struct UnwrapSignal<imj2::Tag, CONTAINER> {
            using T = typename CONTAINER::value_type;
            static auto run(CONTAINER const & container, int N) {
                assert(container.end() == container.begin() + N);
                return complexify<T>(container.begin(), container.begin() + N);
            }
        };
    } // NS fft

    namespace imj2 {
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
    } // NS imj2
} // NS imajuscule
