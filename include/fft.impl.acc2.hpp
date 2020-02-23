#ifndef IMJ_USE_SLOW_FFT

namespace imajuscule {
    // implementation of Accelerate vDSP fft

    namespace accelerate2 {
        struct Tag {};

        static inline unsigned int interlacedIndex(unsigned int b,
                                            unsigned int halfN) {
            assert(halfN);
            unsigned int halfB = b/2;
            unsigned int remainder = b - 2*halfB;
            assert(remainder == 0 || remainder == 1);
            if(remainder) {
                return halfN+halfB;
            }
            else {
                return halfB;
            }
        }
    }

    namespace fft {

        /*
         * Space complexity, for forward fft of real input of size N:
         *
         * input  :  N
         * output :  N
         */

        template<typename T>
        struct RealSignal_<accelerate2::Tag, T> {
            using This = RealSignal_<accelerate2::Tag, T>;
            using type = a64::vector<T>;
            using outputType = type;
            using iter = typename type::iterator;
            using const_iter = typename type::const_iterator;
            using value_type = typename type::value_type;

            static type make(type reals) {
                return std::move(reals);
            }

            static T get_signal(T r) {
                return r;
            }
            
            static void add_assign(value_type * __restrict res,
                                   value_type const * const __restrict const_add,
                                   unsigned int const addSz,
                                   unsigned int const start,
                                   unsigned int const N) {
                using accelerate2::interlacedIndex;
                unsigned int const halfAddSz = [addSz]() {
                    if constexpr (overlapMode == Overlap::Add) {
                        assert(addSz >= 2);
                        return addSz/2;
                    }
                    else {
                        assert(addSz >= 1);
                        return addSz;
                    }
                }();

                for(unsigned int i=start, end=start+N; i!=end; ++i, ++res) {
                    *res += const_add[interlacedIndex(i, halfAddSz)];
                }
            }
            
            
            template<typename TDest>
            static void add_assign_output(TDest * __restrict res,
                                          value_type const * const __restrict const_add,
                                          int N) {
                if constexpr(std::is_same_v<value_type, TDest>) {
                    accelerate::API<T>::f_vadd(res, 1,
                                               const_add, 1,
                                               res, 1,
                                               N);
                }
                else {
                    for(int i=0; i!= N; ++i) {
                        res[i] += const_add[i];
                    }
                }
            }
            
            static void copy(value_type * __restrict dest,
                             value_type const * const __restrict from,
                             unsigned int N) {
                // these 2 are equivalent:
                /*accelerate::API<T>::f_vcpy(N,
                                           &*from, 1,
                                           &*dest, 1);
                */
                accelerate::API<T>::f_mmov(from, dest, 1, N, 1, 1);
            }
            
            template<typename TDest>
            static void copyOutputToOutput(TDest * __restrict dest,
                                           value_type const * const __restrict src,
                                           unsigned int N) {
                if constexpr(std::is_same_v<TDest, value_type>) {
                    copy(dest, src, N);
                }
                else {
                    for(unsigned int i=0; i!= N; ++i) {
                        dest[i] = src[i];
                    }
                }
            }
            
            template<typename TSource>
            static void copyFromInput(value_type * __restrict dest,
                                      TSource const * const __restrict src,
                                      unsigned int N) {
                if constexpr(std::is_same_v<value_type, TSource>) {
                    copy(dest, src, N);
                }
                else {
                    for(unsigned int i=0; i!= N; ++i) {
                        dest[i] = src[i];
                    }
                }
            }
            
            static void copyToOutput(T * __restrict dest,
                                     value_type const * __restrict from,
                                     unsigned int const fromSz,
                                     unsigned int const start,
                                     unsigned int const N) {
                
                unsigned int const halfFromSz = [fromSz]() {
                    if constexpr (overlapMode == Overlap::Add) {
                        assert(fromSz >= 2);
                        return fromSz/2;
                    }
                    else {
                        assert(fromSz >= 1);
                        return fromSz;
                    }
                }();
                
                using accelerate2::interlacedIndex;
                
                for(unsigned int i=start, end=start+N; i!=end; ++i, ++dest) {
                    unsigned int b = interlacedIndex(i, halfFromSz);
                    assert(i < 2*halfFromSz);
                    assert(b < 2*halfFromSz);
                    //assert(b == interlacedIndex(i + 2*halfFromSz, halfFromSz));

                    *dest = from[b];
                }
            }
            
            static void zero_n_raw(T * p, unsigned int n) {
                T zero{};

                accelerate::API<T>::f_vfill(&zero,
                                            p, 1, n);
            }
            static constexpr auto zero_n_raw_output = zero_n_raw;
            
            static void zero_n(type & v, unsigned int n) {
                zero_n_raw(&v[0], n);
            }
            static void zero(type & v) {
                zero_n(v, v.size());
            }
            
            static void dotpr(T const * const a, T const * const b, T * res, unsigned int n) {
                accelerate::API<T>::f_dotpr(a, 1, b, 1, res, n);
            }
        };

        template<typename T, template<typename> typename Allocator>
        struct RealFBins_<accelerate2::Tag, T, Allocator> {
            using This = RealFBins_<accelerate2::Tag, T, Allocator>;
            using Tag = accelerate2::Tag;
            using type = RealFBinsImpl<T, Allocator>;

            // this is slow, it is used for tests only
            static type make(std::vector<complex<T>> const & cplx) {
                // 'wrap' signal
                type res(cplx.size());
                auto split = res.get_hybrid_split();
                *split.realp = cplx[0].real();
                *split.imagp = cplx[cplx.size()/2].real();

                for(int i=1; i<cplx.size()/2; ++i) {
                    split.realp[i] = cplx[i].real();
                    split.imagp[i] = cplx[i].imag();
                }
                return std::move(res);
            }
            
            static void scale(type & v, T const factor) {
                // v *= scalar

                auto V = v.get_hybrid_split();

                accelerate::API<T>::f_vsmul(V.realp, 1,
                                            &factor,
                                            V.realp, 1,
                                            v.size());
            }
            
            static void mult_assign(type & v, type const & const_w) {
                // v *= w

                auto & w = const_cast<type &>(const_w);

                auto V = v.get_hybrid_split();
                auto W = w.get_hybrid_split();

              // handle the first element separately
                *V.realp *= *W.realp;
                *V.imagp *= *W.imagp;

                advance(V);
                advance(W);

                auto const sz = v.vector_size();
                Assert(sz >= 0);
                if(likely(sz > 0)) {
                    accelerate::API<T>::f_zvmul(&V, 1,
                                                &W, 1,
                                                &V, 1,
                                                sz - 1,
                                                1);
                }
            }
            
            static void zero(type & v) {
                T zero{};

                accelerate::SplitComplex<T> sc{
                    &zero,
                    &zero
                };

                auto V = v.get_hybrid_split();
                accelerate::API<T>::f_zvfill(&sc,
                                             &V,
                                             1,
                                             v.vector_size());
            }
            
           static void multiply(T * res,
                                T const * const_m1,
                                T const * const_m2,
                                int N)
            {
                // res = m1 * m2
                
                Assert(N >= 1);
                {
                    res[0] = const_m1[0] * const_m2[0];
                    res[N] = const_m1[N] * const_m2[N];
                }

                if(likely(N > 1)) {
                    using namespace accelerate;

                    SplitComplex<T> Rs {
                        &res[1],
                        &res[N+1]
                    };
                    SplitComplex<T> M1 {
                        const_cast<T *>(&const_m1[1]),
                        const_cast<T *>(&const_m1[N+1])
                    };
                    SplitComplex<T> M2 {
                        const_cast<T *>(&const_m2[1]),
                        const_cast<T *>(&const_m2[N+1])
                    };

                    constexpr bool conjugation = false;
                    accelerate::API<T>::f_zvmul(&M1, 1,
                                                &M2, 1,
                                                &Rs, 1,
                                                N - 1,
                                                conjugation);
                }
            }

            static void multiply_add(T * __restrict accum,
                                     T const * __restrict const_m1,
                                     T const * __restrict const_m2,
                                     int N) {
                // accum += m1 * m2

                Assert(N >= 1);
                {
                    accum[0] += const_m1[0] * const_m2[0];
                    accum[N] += const_m1[N] * const_m2[N];
                }

                if(likely(N > 1)) {
                    using namespace accelerate;

                    SplitComplex<T> Ac {
                        &accum[1],
                        &accum[N+1]
                    };
                    SplitComplex<T> M1 {
                        const_cast<T *>(&const_m1[1]),
                        const_cast<T *>(&const_m1[N+1])
                    };
                    SplitComplex<T> M2 {
                        const_cast<T *>(&const_m2[1]),
                        const_cast<T *>(&const_m2[N+1])
                    };

                    API<T>::f_zvma(&M1, 1,
                                   &M2, 1,
                                   &Ac, 1,
                                   &Ac, 1,
                                   N - 1);
                }
            }

            static std::pair<int, T> getMaxSquaredAmplitude(type const & const_v) {
                typename type::value_type Max = 0;

                auto & v = const_cast<type &>(const_v);
                auto V = v.get_hybrid_split();
                auto index = 0;
                auto const sz = v.vector_size();
                Max = std::max(Max, *V.realp * *V.realp);
                {
                    auto M = *V.imagp * *V.imagp;
                    if(M > Max) {
                        index = sz;
                        Max = M;
                    }
                }

                if(sz > 0) {
                    int const n = sz - 1;
                    for(int i=0; i<n; ++i) {
                        advance(V);
                        auto M = (*V.realp * *V.realp) + (*V.imagp * *V.imagp);
                        if(M > Max) {
                            index = i+1;
                            Max = M;
                        }
                    }
                }

                auto div = static_cast<T>(const_v.size()) * Algo_<Tag,T>::scale;

                return {index, Max/(div * div)};
            }
            
            static void setFirstReal(type & v, T value) {
                if(v.empty()) {
                    throw std::logic_error("setFirstReal on empty");
                }
                v.data()[0] = value;
            }
            
            static T getFirstReal(type const & const_v) {
                auto & v = const_cast<type &>(const_v);
                if(v.empty()) {
                    throw std::logic_error("getFirstReal on empty");
                }
                return v.data()[0];
            }

        };

        template<typename T>
        struct Context_<accelerate2::Tag, T> {
            using type = accelerate::FFTSetup_<T>;

            static auto create(int size) {
                return accelerate::API<T>::f_create_fftsetup(power_of_two_exponent(size), kFFTRadix2);
            }
            static constexpr auto destroy = accelerate::API<T>::f_destroy_fftsetup;
        };

        template<typename T>
        struct Algo_<accelerate2::Tag, T> {
            static constexpr bool inplace_dft = true;

            // it's not clear what I should use :
            // on ios it seems to be a little faster with a tmp buffer,
            // but it's the opposite on osx

            //static constexpr auto ffttype = FFTType::WithTmpBuffer;
            static constexpr auto ffttype = FFTType::Normal;

            using FPT = T;
            using RealInput  = typename RealSignal_ <accelerate2::Tag, T>::type;
            using RealValue  = typename RealInput::value_type;

            template<template<typename> typename Allocator>
            using RealFBins  = typename RealFBins_<accelerate2::Tag, T, Allocator>::type;
            
            using Context    = typename Context_   <accelerate2::Tag, T>::type;
            using Contexts = fft::Contexts_<accelerate2::Tag, T>;
            using Tr = NumTraits<T>;

            // https://developer.apple.com/library/content/documentation/Performance/Conceptual/vDSP_Programming_Guide/UsingFourierTransforms/UsingFourierTransforms.html#//apple_ref/doc/uid/TP40005147-CH202-16195
            static constexpr auto scale = Tr::two();

            Algo_() = default;
            Algo_(Context c) : context(c) {}

            void setContext(Context c) {
                context = c;
            }

          void forward(typename RealInput::const_iterator inputBegin,
                         T * output,
                         unsigned int N) const
            {
                using namespace accelerate;
                SplitComplex<T> Output {
                    output,
                    output + N/2
                };

                constexpr auto inputStride = 1;
                API<T>::f_ctoz(reinterpret_cast<Complex<T> const *>(inputBegin.base()),
                               inputStride * 2,
                               &Output,
                               1,
                               N/2);

                if constexpr (ffttype == FFTType::WithTmpBuffer) {
                    auto & buffer = getFFTTmp();
                    buffer.reserve(N * sizeof(T));

                    SplitComplex<T> buf {
                        reinterpret_cast<T*>(&buffer[0]),
                        reinterpret_cast<T*>(&buffer[0]) + N / 2
                    };

                    API<T>::f_fft_zript(context,
                                        &Output,
                                        1,
                                        &buf,
                                        power_of_two_exponent(N),
                                        FFT_FORWARD);
                }
                else {
                    API<T>::f_fft_zrip(context,
                                       &Output,
                                       1,
                                       power_of_two_exponent(N),
                                       FFT_FORWARD);
                }
            }

            void inverse(T const * const_output,
                         unsigned int N) const
            {
                using namespace accelerate;

                SplitComplex<T> Output {
                    const_cast<T*>(const_output),
                    const_cast<T*>(const_output + N/2)
                };

                if constexpr (ffttype == FFTType::WithTmpBuffer) {
                    auto & buffer = getFFTTmp();
                    buffer.reserve(N * sizeof(T));

                    SplitComplex<T> buf {
                        reinterpret_cast<T*>(&buffer[0]),
                        reinterpret_cast<T*>(&buffer[0]) + N / 2
                    };

                    API<T>::f_fft_zript(context,
                                        &Output,
                                        1,
                                        &buf,
                                        power_of_two_exponent(N),
                                        FFT_INVERSE);
                }
                else {
                    API<T>::f_fft_zrip(context,
                                       &Output,
                                       1,
                                       power_of_two_exponent(N),
                                       FFT_INVERSE);
                }
            }
                        
            static T extractRealOutput(T const * x,
                                       unsigned int b,
                                       unsigned int N)
            {
                using accelerate2::interlacedIndex;
                return x[interlacedIndex(b,
                                         N/2)];
            }

            Context context;
        };

        namespace slow_debug {

            template<typename CONTAINER>
            struct UnwrapFrequenciesRealFBins<accelerate2::Tag, CONTAINER> {
                using T = typename CONTAINER::value_type;
                static auto run(CONTAINER const & const_container, int N) {

                    auto observed = const_cast<CONTAINER &>(const_container).get_hybrid_split();

                    std::vector<complex<T>> res(N, {0,0});
                    res[0] = {
                        observed.realp[0],
                        0
                    };
                    for(int i=1; i<N/2; ++i) {
                        res[i] = {
                            observed.realp[i],
                            observed.imagp[i]
                        };
                    }
                    res[N/2] = {
                        observed.imagp[0],
                        0
                    };
                    const auto pivot = N/2;
                    for(int i=1; i<N/2; ++i) {
                        res[pivot + i] = {
                            +res[pivot - i].real(),
                            -res[pivot - i].imag()
                        };
                    }
                    return std::move(res);
                }
            };

            template<typename CONTAINER>
            struct UnwrapSignal<accelerate2::Tag, CONTAINER> {
                using T = typename CONTAINER::value_type;
                static auto run(CONTAINER const & container, int N) {
                    assert(container.end() == container.begin() + N);
                    return complexify<T>(container.begin(), container.begin() + N);
                }
            };
        } // NS slow_debug
    }// NS fft

    namespace accelerate2 {
        namespace fft {
            using namespace imajuscule::fft;

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
    }// NS accelerate2
}// NS imajuscule

#endif
