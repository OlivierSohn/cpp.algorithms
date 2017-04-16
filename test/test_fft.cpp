
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
        
        template<typename Tag, typename T>
        void testFFT() {
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
            
            Algo fft_algo(setup.get());
            
            fft_algo.run(input, output, N);
            
            // testing against ref found here : https://rosettacode.org/wiki/Fast_Fourier_transform
            
            constexpr auto ffteps = getFFTEpsilon<T>(N);
            
            // we use ffteps when we have an exact value to compare with:
            
            auto res = slow_debug::unwrap<Tag>(output, N);
            
            constexpr auto scale = Algo::scale;
            
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
        
        template<typename Tag>
        void testFFT() {
            testFFT<Tag, float>();
            testFFT<Tag, double>();
        }
    }
}

TEST(FFT, correctness) {
    using namespace imajuscule;
    using namespace imajuscule::testfft;

    testFFT<imj::Tag>();

#if __APPLE__

    testFFT<accelerate::Tag>();

#endif // __APPLE__
    
}
