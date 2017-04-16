
/*
 * Comparative analysis of space complexity, for forward fft of real input of size N.
 *
 *   Summary : AccelerateFFT has 2x better cache locality
 *
 * - MyFFT :
 *
 * roots :   N
 * input : 2*N
 * result: 2*N
 *
 * - AccelerateFFT:
 *
 * setup :   N (?)
 * input :   N
 * result:   N
 *
 */

namespace imajuscule {
    namespace testfft {
        
        template<typename T>
        constexpr auto getFFTEpsilon(int N) {
            return power_of_two_exponent(N) * std::numeric_limits<T>::epsilon(); // worst case error propagation is O(log N)
        }
        
        template<typename T>
        void createRealInput(T & input) {
            for(int i=0; i<4; ++i) {
                input.emplace_back(1);
            }
            for(int i=0; i<4; ++i) {
                input.emplace_back(0);
            }
        }
        
        template<typename Tag, typename T, typename Output>
        void verifyFrequencies(Output output, int N) {
            using namespace imajuscule::fft;
            
            // testing against ref found here : https://rosettacode.org/wiki/Fast_Fourier_transform
            
            const auto ffteps = getFFTEpsilon<T>(N);
            
            // we use ffteps when we have an exact value to compare with:
            
            auto const res = slow_debug::unwrap_frequencies<Tag>(output, N);
            
            constexpr auto scale = Algo_<Tag, T>::scale;
            
            ASSERT_NEAR(scale * 0, res[0].imag(), ffteps);
            ASSERT_NEAR(scale * 4, res[0].real(), ffteps);
            ASSERT_NEAR(scale * 0, res[2].imag(), ffteps);
            ASSERT_NEAR(scale * 0, res[2].real(), ffteps);
            ASSERT_NEAR(scale * 0, res[4].imag(), ffteps);
            ASSERT_NEAR(scale * 0, res[4].real(), ffteps);
            ASSERT_NEAR(scale * 0, res[6].imag(), ffteps);
            ASSERT_NEAR(scale * 0, res[6].real(), ffteps);
            
            ASSERT_NEAR(scale * 1, res[1].real(), ffteps);
            ASSERT_NEAR(scale * -2.41421, res[1].imag(), 1e-4);
            ASSERT_NEAR(scale * 1, res[3].real(), ffteps);
            ASSERT_NEAR(scale * -0.414214, res[3].imag(), 1e-4);
            ASSERT_NEAR(scale * 1, res[5].real(), ffteps);
            ASSERT_NEAR(scale * 0.414214, res[5].imag(), 1e-4);
            ASSERT_NEAR(scale * 1, res[7].real(), ffteps);
            ASSERT_NEAR(scale * 2.41421, res[7].imag(), 1e-4);
        }
        
        template<typename Tag, typename T, typename Input>
        void verifySignal(Input signal, int N) {
            using namespace imajuscule::fft;
            
            const auto ffteps = getFFTEpsilon<T>(N);
            
            // we use ffteps when we have an exact value to compare with:
                        
            auto const res = slow_debug::unwrap_signal<Tag>(signal, N);
            
            constexpr auto scale = 1;
            
            ASSERT_NEAR(scale * 1, res[0].real(), ffteps);
            ASSERT_NEAR(scale * 0, res[0].imag(), ffteps);
            ASSERT_NEAR(scale * 1, res[1].real(), ffteps);
            ASSERT_NEAR(scale * 0, res[1].imag(), ffteps);
            ASSERT_NEAR(scale * 1, res[2].real(), ffteps);
            ASSERT_NEAR(scale * 0, res[2].imag(), ffteps);
            ASSERT_NEAR(scale * 1, res[3].real(), ffteps);
            ASSERT_NEAR(scale * 0, res[3].imag(), ffteps);
            ASSERT_NEAR(scale * 0, res[4].real(), ffteps);
            ASSERT_NEAR(scale * 0, res[4].imag(), ffteps);
            ASSERT_NEAR(scale * 0, res[5].real(), ffteps);
            ASSERT_NEAR(scale * 0, res[5].imag(), ffteps);
            ASSERT_NEAR(scale * 0, res[6].real(), ffteps);
            ASSERT_NEAR(scale * 0, res[6].imag(), ffteps);
            ASSERT_NEAR(scale * 0, res[7].real(), ffteps);
            ASSERT_NEAR(scale * 0, res[7].imag(), ffteps);
        }

        template<typename Tag, typename T>
        void testForwardFFT() {
            using namespace imajuscule;
            using namespace imajuscule::fft;
            using namespace imajuscule::testfft;

            using RealInput = typename RealInput_<Tag, T>::type;
            using RealOutput = typename RealOutput_<Tag, T>::type;
            using Context = typename Context_<Tag, T>::type;
            using ScopedContext = ScopedContext_<Tag, T>;
            using Algo = Algo_<Tag, T>;

            constexpr auto N = 8;
            ScopedContext setup(N);
            
            RealInput input;
            RealOutput output(N);
            
            input.reserve(N);
            
            createRealInput(input);

            { verifySignal<Tag, T>(input, N); }

            Algo fft_algo(setup.get());
            
            fft_algo.forward(input, output, N);
            
            { verifyFrequencies<Tag, T>(output, N); }
        }
        
        template<typename Tag>
        void testForwardFFT() {
            testForwardFFT<Tag, float>();
            testForwardFFT<Tag, double>();
        }
        
        template<typename Tag, typename T>
        void testInverseFFT() {
            using namespace imajuscule;
            using namespace imajuscule::fft;
            using namespace imajuscule::testfft;
            
            using RealInput = typename RealInput_<Tag, T>::type;
            using RealOutput = typename RealOutput_<Tag, T>::type;
            using Context = typename Context_<Tag, T>::type;
            using ScopedContext = ScopedContext_<Tag, T>;
            using Algo = Algo_<Tag, T>;
            
            constexpr auto N = 8;
            ScopedContext setup(N);
            
            RealInput input, reconstructed_input(N);
            RealOutput output(N);
            
            input.reserve(N);
            
            createRealInput(input);
            
            { verifySignal<Tag, T>(input, N); }
            
            Algo fft_algo(setup.get());
            
            fft_algo.forward(input, output, N);
            
            { verifyFrequencies<Tag, T>(output, N); }
            
            fft_algo.inverse(output, reconstructed_input, N);
            
            for(auto & v : reconstructed_input) {
                v *= 1/(Algo::scale * static_cast<T>(N));
            }
          
            { verifySignal<Tag, T>(reconstructed_input, N); }
            
        }
        
        template<typename Tag>
        void testInverseFFT() {
            testInverseFFT<Tag, float>();
            testInverseFFT<Tag, double>();
        }
    }
}

TEST(FFT, forward_correctness) {
    using namespace imajuscule;
    using namespace imajuscule::testfft;
    
    testForwardFFT<imj::Tag>();
    
#if __APPLE__
    
    testForwardFFT<accelerate::Tag>();
    
#endif // __APPLE__
}

TEST(FFT, inverse_correctness) {
    using namespace imajuscule;
    using namespace imajuscule::testfft;
    
    testInverseFFT<imj::Tag>();
    
#if __APPLE__
    
    testInverseFFT<accelerate::Tag>();
    
#endif // __APPLE__
}
