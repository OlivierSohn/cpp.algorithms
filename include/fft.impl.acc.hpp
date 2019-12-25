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
            
            static double cost_add_scalar_multiply(int64_t sz) {
                static std::unordered_map<int, double> stats;
                static std::mutex mut;
                std::lock_guard<std::mutex> l(mut);
                
                auto it = stats.find(sz);
                if(it != stats.end()) {
                    return it->second;
                }
                double t = nocache_cost_add_scalar_multiply(sz);
                stats[sz] = t;
                return t;
            }
            
            static double nocache_cost_add_scalar_multiply(int64_t sz) {
                std::array<type, 3> a;
                for(auto & v:a) {
                    v.resize(sz, 1.);
                }
                int64_t ntests = std::max(static_cast<int64_t>(1),
                                          10000 / sz);
                using namespace profiling;
                auto duration = measure_thread_cpu_one([&a, sz, ntests](){
                    std::array<int, 3> idx {0,1,2};
                    
                    for(int i=0; i<ntests; i++) {
                        This::add_scalar_multiply(a[idx[0]].begin(),
                                                  a[idx[1]].begin(),
                                                  a[idx[2]].begin(),
                                                  0.7,
                                                  sz);
                        std::rotate(idx.begin(),
                                    idx.begin()+1,
                                    idx.end());
                    }
                });
                T avg = std::abs(std::accumulate(a[2].begin(), a[2].end(), T{})) / static_cast<T>(a[2].size());
                return duration.count() / (1000000. * static_cast<double>(ntests) + std::min(0.000001, avg));
            }
            
            static void copy(iter dest, const_iter from, int N) {
                // these 2 are equivalent:
                /*accelerate::API<T>::f_vcpy(N,
                                           &*from, 1,
                                           &*dest, 1);
                */
                accelerate::API<T>::f_mmov(&*from, &*dest, 1, N, 1, 1);
            }
            
            static double cost_copy(int64_t sz) {
                static std::unordered_map<int, double> stats;
                static std::mutex mut;
                std::lock_guard<std::mutex> l(mut);
                
                auto it = stats.find(sz);
                if(it != stats.end()) {
                    return it->second;
                }
                double t = nocache_cost_copy(sz);
                stats[sz] = t;
                return t;
            }
            
            static auto nocache_cost_copy(int64_t sz) {
                std::array<type, 2> a;
                {
                    int i=0;
                    for(auto & v:a) {
                        ++i;
                        v.resize(sz, static_cast<T>(i));
                    }
                }
                int64_t ntests = std::max(static_cast<int64_t>(1),
                                          10000 / sz);
                using namespace profiling;
                auto duration = measure_thread_cpu_one([&a, sz, ntests](){
                    std::array<int, 3> idx {0,1};
                    for(int i=0; i<ntests; i++) {
                        This::copy(a[idx[0]].begin(),
                                   a[idx[1]].begin(),
                                   sz);
                        std::rotate(idx.begin(),
                                    idx.begin()+1,
                                    idx.end());
                    }
                });
                T avg = std::abs(std::accumulate(a[1].begin(), a[1].end(), T{})) / static_cast<T>(a[1].size());
                return duration.count() / (1000000. * static_cast<double>(ntests) + std::min(0.000001, avg));
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
          auto empty() const { return buffer.empty(); }

            auto vector_size() const {
                return buffer.size() / 2;
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
            a64::vector<T> buffer;
        };

        template<typename ComplexSplit>
        void advance(ComplexSplit & cs) {
            ++cs.realp;
            ++cs.imagp;
        }

        template<typename T>
        struct RealFBins_<accelerate::Tag, T> {
            using This = RealFBins_<accelerate::Tag, T>;
            using Tag = accelerate::Tag;
            using type = RealFBinsImpl<T>;

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
            
            static double cost_mult_assign(int64_t sz) {
                static std::unordered_map<int, double> stats;
                static std::mutex mut;
                std::lock_guard<std::mutex> l(mut);
                
                auto it = stats.find(sz);
                if(it != stats.end()) {
                    return it->second;
                }
                double t = nocache_cost_mult_assign(sz);
                stats[sz] = t;
                return t;
            }
            
            static auto nocache_cost_mult_assign(int64_t sz) {
                std::array<type, 2> a;
                int i=0;
                for(auto & v:a) {
                    ++i;
                    v.resize(sz);
                    *v.get_hybrid_split().realp = i;
                }
                int64_t ntests = std::max(static_cast<int64_t>(1),
                                          10000 / sz);
                using namespace profiling;
                auto duration = measure_thread_cpu_one([&a, sz, ntests](){
                    for(int i=0; i<ntests; ++i) {
                        This::mult_assign(a[0],
                                          a[1]);
                    }
                });
                T sum = 0;
                auto s = a[0].get_hybrid_split();
                for(int i=0; i<a[0].vector_size(); ++i) {
                    sum += s.realp[i];
                    sum += s.imagp[i];
                }
                sum /= a[0].size();
                sum = std::abs(sum);
                return duration.count() / (1000000. * static_cast<double>(ntests) + std::min(0.000001, sum));
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

                auto const sz = accum.vector_size();
                if(likely(sz > 0)) {
                    accelerate::API<T>::f_zvma(&M1, 1,
                                               &M2, 1,
                                               &Accum, 1,
                                               &Accum, 1, sz - 1);
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

            using RealInput  = typename RealSignal_ <accelerate::Tag, T>::type;
            using RealFBins  = typename RealFBins_<accelerate::Tag, T>::type;
            using Context    = typename Context_   <accelerate::Tag, T>::type;
            using Contexts = fft::Contexts_<accelerate::Tag, T>;
            using Tr = NumTraits<T>;

            // scaling factor of 2 :
            // https://developer.apple.com/library/content/documentation/Performance/Conceptual/vDSP_Programming_Guide/UsingFourierTransforms/UsingFourierTransforms.html#//apple_ref/doc/uid/TP40005147-CH202-16195
            static constexpr auto scale = Tr::two();

            Algo_() = default;
            Algo_(Context c) : context(c) {}

            void setContext(Context c) {
                context = c;
            }

          void forward(typename RealInput::const_iterator inputBegin,
                         RealFBins & output,
                         unsigned int N) const
            {
                using namespace accelerate;
                auto Output = output.get_hybrid_split();

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

            void inverse(RealFBins const & const_output,
                         RealInput & input,
                         unsigned int N) const
            {
                using namespace accelerate;

                auto Output = const_cast<RealFBins &>(const_output).get_hybrid_split();

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
                               reinterpret_cast<Complex<T> *>(&input[0]),
                               inputStride * 2,
                               N/2);

            }
            
            static double cost_fft_inverse(int64_t sz) {
                static std::unordered_map<int, double> stats;
                static std::mutex mut;
                std::lock_guard<std::mutex> l(mut);
                
                auto it = stats.find(sz);
                if(it != stats.end()) {
                    return it->second;
                }
                double t = nocache_cost_fft_inverse(sz);
                stats[sz] = t;
                return t;
            }
            
            
            static auto nocache_cost_fft_inverse(int64_t sz) {
                RealFBins f;
                RealInput input;
                f.resize(sz);
                input.resize(sz);
                Algo_<accelerate::Tag, T> a(Contexts::getInstance().getBySize(sz));

                int64_t ntests = std::max(static_cast<int64_t>(1),
                                          10000 / sz);

                using namespace profiling;
                T sum = 0;
                auto duration = measure_thread_cpu_one([&sum, &a, &f, &input, sz, ntests](){
                    for(int i=0; i<ntests; ++i) {
                        *f.get_hybrid_split().realp = static_cast<T>(i);
                        a.inverse(f, input, sz);
                        sum += input[0];
                    }
                });
                sum = std::abs(sum);
                return duration.count() / (1000000. * static_cast<double>(ntests) + std::min(0.000001, sum));
            }
                        
            static double cost_fft_forward(int64_t sz) {
                static std::unordered_map<int, double> stats;
                static std::mutex mut;
                std::lock_guard<std::mutex> l(mut);
                
                auto it = stats.find(sz);
                if(it != stats.end()) {
                    return it->second;
                }
                double t = nocache_cost_fft_forward(sz);
                stats[sz] = t;
                return t;
            }
            
            static auto nocache_cost_fft_forward(int64_t sz) {
                RealFBins f;
                RealInput input;
                f.resize(sz);
                input.resize(sz);
                Algo_<accelerate::Tag, T> a(Contexts::getInstance().getBySize(sz));
                
                int64_t ntests = std::max(static_cast<int64_t>(1),
                                          10000 / sz);

                using namespace profiling;
                T sum = 0;
                auto duration = measure_thread_cpu_one([&sum, &a, &f, &input, sz, ntests](){
                    for(int i=0; i<ntests; ++i) {
                        input[0] = static_cast<T>(i);
                        a.forward(input.begin(), f, sz);
                        sum += *f.get_hybrid_split().realp;
                    }
                });
                sum = std::abs(sum);
                return duration.count() / (1000000. * static_cast<double>(ntests) + std::min(0.000001, sum));
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

#endif
