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
         * input :   N
         */

        template<typename T>
        struct RealInput_<accelerate::Tag, T> {
            using type = std::vector<T>;
        };
        

        template<typename T>
        struct RealOutputImpl {
            using value_type = T;
            
            RealOutputImpl(int size) : observed_mem(size) {
                observed = { observed_mem.data(), observed_mem.data() + size/2 };
            }
            
            accelerate::SplitComplex<T> observed;
            private:
            std::vector<T> observed_mem;
        };
        
        template<typename T>
        struct RealOutput_<accelerate::Tag, T> {
            using Tag = accelerate::Tag;
            using type = RealOutputImpl<T>;
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
            using RealInput  = typename RealInput_ <accelerate::Tag, T>::type;
            using RealOutput = typename RealOutput_<accelerate::Tag, T>::type;
            using Context    = typename Context_   <accelerate::Tag, T>::type;
            
            Algo_(Context c) : ctxt(c) {}
            
            void run(RealInput const & input,
                     RealOutput & output,
                     unsigned int N) const
            {
                using namespace accelerate;
                auto & observed = output.observed;
                {
                    constexpr auto inputStride = 1;
                    ctoz<T>(reinterpret_cast<Complex<T> const *>(input.data()),
                            inputStride * 2,
                            &observed,
                            1,
                            N/2);
                    
                    fft_zrip<T>(ctxt,
                                &observed,
                                1,
                                power_of_two_exponent(N),
                                FFT_FORWARD);
                }
            }
            Context ctxt;
        };
        
        namespace slow_debug {
            template<typename CONTAINER>
            struct Unwrap<accelerate::Tag, CONTAINER> {
                using T = typename CONTAINER::value_type;
                static auto run(CONTAINER const & container, int N) {
                    
                    auto & observed = container.observed;
                    
                    FFTVec<T> res(N, {0,0});
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
        }
    }
    
    namespace accelerate {
        namespace fft {
            using namespace imajuscule::fft;
            
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
        }
    }
}

