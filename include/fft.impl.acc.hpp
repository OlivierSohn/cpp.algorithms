#ifndef IMJ_USE_SLOW_FFT

namespace imajuscule {
    // implementation of Accelerate vDSP fft

    namespace accelerate {
        struct Tag {};
    }

    namespace fft {

        /*
         * Space complexity, for forward fft of real input of size N:
         *
         * input  :  N
         * output :  N
         */

        template<typename T>
        struct RealSignal_<accelerate::Tag, T> {
            using This = RealSignal_<accelerate::Tag, T>;
            using type = a64::vector<T>;
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
                                   int N) {
                // res += add

                accelerate::API<T>::f_vadd(res, 1,
                                           const_add, 1,
                                           res, 1,
                                           N);
            }
            
            static void add_scalar_multiply(iter res,
                                            const_iter const_add1,
                                            const_iter const_add2, T const m, int N) {
                // res = m * (add1 + add2)

                accelerate::API<T>::f_vasm(&*const_add1, 1,
                                           &*const_add2, 1,
                                           &m,
                                           &*res, 1,
                                           N);

            }
            
            static void copy(value_type * __restrict dest,
                             value_type const * const __restrict from,
                             int N) {
                // these 2 are equivalent:
                /*accelerate::API<T>::f_vcpy(N,
                                           &*from, 1,
                                           &*dest, 1);
                */
                accelerate::API<T>::f_mmov(from, dest, 1, N, 1, 1);
            }
            
            static void zero_n_raw(T * p, int n) {
                T zero{};

                accelerate::API<T>::f_vfill(&zero,
                                            p, 1, n);
            }            
            static void zero_n(type & v, int n) {
                zero_n_raw(&v[0], n);
            }
            static void zero(type & v) {
                zero_n(v, v.size());
            }
            
            static void dotpr(T const * const a, T const * const b, T * res, int n) {
                accelerate::API<T>::f_dotpr(a, 1, b, 1, res, n);
            }
        };


        /*
         Represents the first half of the spectrum (the second half is the conjugate)
         and the nyquist frequency (real) is encoded in the 0th index imag.
         cf. packing here : https://developer.apple.com/library/content/documentation/Performance/Conceptual/vDSP_Programming_Guide/UsingFourierTransforms/UsingFourierTransforms.html
         */
        template<typename T, template<typename> typename A>
        struct RealFBinsImpl {
            using SC = accelerate::SplitComplex<T>;

            using value_type = T;
            using allocator_type = A<T>;

            RealFBinsImpl() = default;

            RealFBinsImpl(int size) : buffer(size) {
            }

            void reserve(int sz) {
                buffer.reserve(sz);
            }
            void resize(size_t sz) {
                buffer.resize(sz);
            }
            void clear() {
                buffer.clear();
            }
            
            auto size() const { return buffer.size(); }
            auto empty() const { return buffer.empty(); }
            
            auto vector_size() const {
                return buffer.size() / 2;
            }
            
            auto data() {
                return buffer.data();
            }
            auto data() const {
                return buffer.data();
            }

            auto get_hybrid_split() {
                Assert(buffer.size() >= 2);
                Assert(0 == buffer.size() % 2);
                return SC {
                    // 64 byte aligned:
                    &buffer[0],
                    // for doubles, and a buffer size of at least 4, this is at least 16 bytes aligned:
                    &buffer[0] + buffer.size()/2
                };
            }

        private:
            std::vector<T, allocator_type> buffer;
        };

        template<typename ComplexSplit>
        void advance(ComplexSplit & cs) {
            ++cs.realp;
            ++cs.imagp;
        }

        template<typename T, template<typename> typename Allocator>
        struct RealFBins_<accelerate::Tag, T, Allocator> {
            using This = RealFBins_<accelerate::Tag, T, Allocator>;
            using Tag = accelerate::Tag;
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

                    accelerate::API<T>::f_zvma(&M1, 1,
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
        struct Context_<accelerate::Tag, T> {
            using type = accelerate::FFTSetup_<T>;

            static auto create(int size) {
                return accelerate::API<T>::f_create_fftsetup(power_of_two_exponent(size), kFFTRadix2);
            }
            static constexpr auto destroy = accelerate::API<T>::f_destroy_fftsetup;
        };

        enum class FFTType {
            Normal,
            WithTmpBuffer
        };

        a64::vector<int8_t> & getFFTTmp();

        template<typename T>
        struct Algo_<accelerate::Tag, T> {

            // it's not clear what I should use :
            // on ios it seems to be a little faster with a tmp buffer,
            // but it's the opposite on osx

            //static constexpr auto ffttype = FFTType::WithTmpBuffer;
            static constexpr auto ffttype = FFTType::Normal;

            using FPT = T;
            using RealInput  = typename RealSignal_ <accelerate::Tag, T>::type;
            using RealValue  = typename RealInput::value_type;

            template<template<typename> typename Allocator>
            using RealFBins  = typename RealFBins_<accelerate::Tag, T, Allocator>::type;
            
            using Context    = typename Context_   <accelerate::Tag, T>::type;
            using Contexts = fft::Contexts_<accelerate::Tag, T>;
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
                         RealValue * input,
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

                constexpr auto inputStride = 1;
                API<T>::f_ztoc(&Output,
                               1,
                               reinterpret_cast<Complex<T> *>(input),
                               inputStride * 2,
                               N/2);

            }
            
            Context context;
        };

        namespace slow_debug {

            template<typename CONTAINER>
            struct UnwrapFrequenciesRealFBins<accelerate::Tag, CONTAINER> {
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
            struct UnwrapSignal<accelerate::Tag, CONTAINER> {
                using T = typename CONTAINER::value_type;
                static auto run(CONTAINER const & container, int N) {
                    assert(container.end() == container.begin() + N);
                    return complexify<T>(container.begin(), container.begin() + N);
                }
            };
        } // NS slow_debug
    }// NS fft

    namespace accelerate {
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
    }// NS accelerate
}// NS imajuscule

#endif
