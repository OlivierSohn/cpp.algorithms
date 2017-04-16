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
    }
}

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

template<typename T>
void testMyFFT() {
    using namespace imajuscule;
    using namespace imajuscule::fft;
    using namespace imajuscule::imj::fft;
    using namespace imajuscule::testfft;

    using Algo = imajuscule::imj::fft::Algo<T>;

    constexpr auto N = 8;
    ScopedContext<T> setup(N);
    
    RealInput<T> input;
    RealOutput<T> output(N);

    input.reserve(N);
    
    createRealInput(input);
    
    Algo fft_algo(setup.get());
    
    fft_algo.run(input, output, N);
        
    // testing against ref found here : https://rosettacode.org/wiki/Fast_Fourier_transform
    
    constexpr auto ffteps = getFFTEpsilon<T>(N);
    
    // we use ffteps when we have an exact value to compare with:
    
    auto const & res = output;
    
    ASSERT_NEAR(0, res[0].imag(), ffteps);
    ASSERT_NEAR(4, res[0].real(), ffteps);
    ASSERT_NEAR(0, res[2].imag(), ffteps);
    ASSERT_NEAR(0, res[2].real(), ffteps);
    ASSERT_NEAR(0, res[4].imag(), ffteps);
    ASSERT_NEAR(0, res[4].real(), ffteps);
    ASSERT_NEAR(0, res[6].imag(), ffteps);
    ASSERT_NEAR(0, res[6].real(), ffteps);
    
    ASSERT_NEAR(1, res[1].real(), ffteps);
    ASSERT_NEAR(-2.41421, res[1].imag(), 1e-4);
    ASSERT_NEAR(1, res[3].real(), ffteps);
    ASSERT_NEAR(-0.414214, res[3].imag(), 1e-4);
    ASSERT_NEAR(1, res[5].real(), ffteps);
    ASSERT_NEAR( 0.414214, res[5].imag(), 1e-4);
    ASSERT_NEAR(1, res[7].real(), ffteps);
    ASSERT_NEAR( 2.41421, res[7].imag(), 1e-4);
}

#if __APPLE__

/*
 * space complexity for real input of size N:
 */
template<typename T>
void testAccelerateFFT() {
    using namespace imajuscule;
    using namespace imajuscule::accelerate;
    using namespace imajuscule::fft;
    using namespace imajuscule::accelerate::fft;

    using Algo = imajuscule::accelerate::fft::Algo<T>;
    using TAG = imajuscule::accelerate::Tag;
    
    constexpr auto N = 8;    
    ScopedContext<T> setup(N);
    
    RealInput<T> input;
    RealOutput<T> output(N);
    
    testfft::createRealInput(input);
    
    Algo fft_algo(setup.get());
    
    fft_algo.run(input, output, N);
    
    // put the result in a complex vector to match the format of the other test
  
    // testing against ref found here : https://rosettacode.org/wiki/Fast_Fourier_transform
    
    constexpr auto ffteps = testfft::getFFTEpsilon<T>(N);
    
    // we use ffteps when we have an exact value to compare with:
    
    // there is a scaling factor of 2 :
    // https://developer.apple.com/library/content/documentation/Performance/Conceptual/vDSP_Programming_Guide/UsingFourierTransforms/UsingFourierTransforms.html#//apple_ref/doc/uid/TP40005147-CH202-16195
    
    constexpr auto scale = 2.f;
    
    auto res = slow_debug::unwrap<TAG>(output, N);
    
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

#endif // __APPLE__

TEST(FFT, simple) {
    testMyFFT<float>();
    testMyFFT<double>();
    testAccelerateFFT<float>();
    testAccelerateFFT<double>();
}

template<typename T>
void testBig() {
    using namespace imajuscule;
    using namespace imajuscule::fft;
    
    constexpr auto N = 4096;
    FFTVec<T> res, input;
    
    res.resize(N);
    input.resize(N);
    auto roots = compute_roots_of_unity<T>(N);
    Algo<T> a(std::move(roots));
    a.run(input.begin(), res.begin(), N, 1 );
}

TEST(FFT, big) {
    testBig<float>();
    testBig<double>();
}

