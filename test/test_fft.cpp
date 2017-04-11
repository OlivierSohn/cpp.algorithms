
template<typename T>
void testFFT() {
    using namespace imajuscule;
    using namespace imajuscule::fft;
    
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
    
    constexpr auto eps = std::numeric_limits<T>::epsilon();
    constexpr auto ffteps = power_of_two_exponent(N) * eps; // worst case error propagation is O(log N)
    
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

TEST(FFT, simple) {
    testFFT<float>();
    testFFT<double>();
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

