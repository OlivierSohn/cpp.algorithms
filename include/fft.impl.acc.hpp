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

        /*
         * Space complexity, for forward fft of real input of size N:
         *
         * input  :  N
         * output :  N
         */

        template<typename T>
        struct RealSignal_<accelerate::Tag, T> {
            using type = cacheline_aligned_allocated::vector<T>;
            using iter = typename type::iterator;
            using const_iter = typename type::const_iterator;
        
            static type make(type reals) {
                return std::move(reals);
            }

            static T get_signal(T r) {
                return r;
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
            
            static void copy(iter dest, const_iter from, int N) {
                // don't know which one is faster .. todo test!
                /*accelerate::API<T>::f_vcpy(N,
                                           &*from, 1,
                                           &*dest, 1);*/
                accelerate::API<T>::f_mmov(&*from, &*dest, 1, N, 1, 1);
            }
        };
        

        /*
         Represents the first half of the spectrum (the second half is the conjugate)
         and the nyquist frequency (real) is encoded in the 0th index imag.
         cf. packing here : https://developer.apple.com/library/content/documentation/Performance/Conceptual/vDSP_Programming_Guide/UsingFourierTransforms/UsingFourierTransforms.html
         */
        template<typename T>
        struct RealFBinsImpl {
            using SC = accelerate::SplitComplex<T>;
            
            using value_type = T;
            
            RealFBinsImpl() = default;
            
            RealFBinsImpl(int size) : buffer(size) {
            }
            
            void resize(size_t sz) {
                buffer.resize(sz);
            }
            
            auto size() const { return buffer.size(); }
            
            auto count_cplx_elements() const { return buffer.size() / 2; }
            
            auto get_hybrid_split() {
                return SC {
                    buffer.data(),
                    buffer.data() + buffer.size()/2
                };
            }

        private:
            cacheline_aligned_allocated::vector<T> buffer;
        };
        
        template<typename ComplexSplit>
        void advance(ComplexSplit & cs) {
            ++cs.realp;
            ++cs.imagp;
        }

        template<typename T>
        struct RealFBins_<accelerate::Tag, T> {
            using Tag = accelerate::Tag;
            using type = RealFBinsImpl<T>;
           
            // this is slow, it is used for tests only
            static type make(std::vector<complex<T>> cplx) {
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

            static void mult_assign(type & v, type const & const_w) {
                // v *= w
                
                auto & w = const_cast<type &>(const_w);
                
                auto V = v.get_hybrid_split();
                auto W = w.get_hybrid_split();

                *V.realp *= *W.realp;
                *V.imagp *= *W.imagp;
                
                advance(V);
                advance(W);
                
                accelerate::API<T>::f_zvmul(&V, 1,
                                            &W, 1,
                                            &V, 1, v.count_cplx_elements()-1, 1);
            }
            
            static void zero(type & v) {
                T zero{};
                
                accelerate::SplitComplex<T> sc{
                    &zero,
                    &zero
                };
                
                auto V = v.get_hybrid_split();
                accelerate::API<T>::f_zvfill(&sc,
                                             &V, 1, v.count_cplx_elements());
            }
            
            static void multiply_add(type & accum, type const & const_m1, type const & const_m2) {
                // accum += m1 * m2
                
                auto & m1 = const_cast<type &>(const_m1);
                auto & m2 = const_cast<type &>(const_m2);
                
                auto Accum = accum.get_hybrid_split();
                auto M1 = m1.get_hybrid_split();
                auto M2 = m2.get_hybrid_split();

                {
                    *Accum.realp += *M1.realp * *M2.realp;
                    *Accum.imagp += *M1.imagp * *M2.imagp;
                }
                
                advance(Accum);
                advance(M1);
                advance(M2);
                
                accelerate::API<T>::f_zvma(&M1, 1,
                                           &M2, 1,
                                           &Accum, 1,
                                           &Accum, 1, accum.count_cplx_elements() - 1);
                
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

        template<typename T>
        struct Algo_<accelerate::Tag, T> {
            using RealInput  = typename RealSignal_ <accelerate::Tag, T>::type;
            using RealFBins  = typename RealFBins_<accelerate::Tag, T>::type;
            using Context    = typename Context_   <accelerate::Tag, T>::type;
            using Tr = NumTraits<T>;

            // scaling factor of 2 :
            // https://developer.apple.com/library/content/documentation/Performance/Conceptual/vDSP_Programming_Guide/UsingFourierTransforms/UsingFourierTransforms.html#//apple_ref/doc/uid/TP40005147-CH202-16195
            static constexpr auto scale = Tr::two();
            
            Algo_() = default;
            Algo_(Context c) : context(c) {}
            
            void setContext(Context c) {
                context = c;
            }
            
            void forward(RealInput const & input,
                         RealFBins & output,
                         unsigned int N) const
            {
                using namespace accelerate;
                auto Output = output.get_hybrid_split();
                
                constexpr auto inputStride = 1;
                ctoz<T>(reinterpret_cast<Complex<T> const *>(input.data()),
                        inputStride * 2,
                        &Output,
                        1,
                        N/2);
                
                fft_zrip<T>(context,
                            &Output,
                            1,
                            power_of_two_exponent(N),
                            FFT_FORWARD);
            }
            
            void inverse(RealFBins const & const_output,
                         RealInput & input,
                         unsigned int N) const
            {
                using namespace accelerate;

                auto Output = const_cast<RealFBins &>(const_output).get_hybrid_split();
                
                fft_zrip<T>(context,
                            &Output,
                            1,
                            power_of_two_exponent(N),
                            FFT_INVERSE);

                constexpr auto inputStride = 1;
                ztoc<T>(&Output,
                        1,
                        reinterpret_cast<Complex<T> *>(input.data()),
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
            
            template<typename T>
            using RealFBins = typename RealFBins_<Tag, T>::type;
            
            template<typename T>
            using Context = typename Context_<Tag, T>::type;
            
            template<typename T>
            using ScopedContext = ScopedContext_<Tag, T>;
            
            template<typename T>
            using Algo = Algo_<Tag, T>;
        } // NS fft
    }// NS accelerate
}// NS imajuscule

