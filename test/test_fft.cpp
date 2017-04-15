namespace imajuscule {
    namespace testfft {
        
        template<typename T>
        constexpr auto getFFTEpsilon(int N) {
            return power_of_two_exponent(N) * std::numeric_limits<T>::epsilon(); // worst case error propagation is O(log N)
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
 * =>
 *
 */

template<typename T>
void testMyFFT() {
    using namespace imajuscule;
    using namespace imajuscule::fft;
    using namespace imajuscule::testfft;
    
    constexpr auto N = 8;
    FFTVec<T> res, input{{
        {1,0},
        {1,0},
        {1,0},
        {1,0},
        {0,0},
        {0,0},
        {0,0},
        {0,0}
    }};
    
    res.resize(N);
    input.resize(N);
    auto roots = compute_roots_of_unity<T>(N);
    Algo<T> a(std::move(roots));
    a.run(input.begin(), res.begin(), N, 1 );
    
    // testing against ref found here : https://rosettacode.org/wiki/Fast_Fourier_transform
    
    constexpr auto ffteps = getFFTEpsilon<T>(N);
    
    // we use ffteps when we have an exact value to compare with:
    
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

/*
 * space complexity for real input of size N:
 */
template<typename T>
void testAccelerateFFT() {
    using namespace imajuscule;
    using namespace imajuscule::accelerate;
    using namespace imajuscule::fft;
    using namespace imajuscule::testfft;
    
    constexpr auto N = 8;
    constexpr auto Log2N = power_of_two_exponent(N);
    
    ScopedFFTSetup<T> setup(Log2N, kFFTRadix2);
    
    constexpr int inputStride = 1;
    std::vector<T> input {{
        1,
        1,
        1,
        1,
        0,
        0,
        0,
        0
    }};
    
    std::vector<T> observed_mem( N );
    
    SplitComplex<T> observed = { observed_mem.data(), observed_mem.data() + N/2 };
    
    ctoz<T>(reinterpret_cast<Complex<T> *>(input.data()),
            inputStride * 2,
            &observed,
            1,
            N/2);
    
    fft_zrip<T>(setup.get(),
                &observed,
                1,
                Log2N,
                FFT_FORWARD);
    
    // put the result in a complex vector to match the format of the other test
    
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
    constexpr auto pivot = N/2;
    for(int i=1; i<N/2; ++i) {
        res[pivot + i] = {
            +res[pivot - i].real(),
            -res[pivot - i].imag()
        };
    }
    // testing against ref found here : https://rosettacode.org/wiki/Fast_Fourier_transform
    
    constexpr auto ffteps = getFFTEpsilon<T>(N);
    
    // we use ffteps when we have an exact value to compare with:
    
    // there is a scaling factor of 2 :
    // https://developer.apple.com/library/content/documentation/Performance/Conceptual/vDSP_Programming_Guide/UsingFourierTransforms/UsingFourierTransforms.html#//apple_ref/doc/uid/TP40005147-CH202-16195
    
    constexpr auto scale = 2.f;
    
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

