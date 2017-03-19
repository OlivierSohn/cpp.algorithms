

TEST(FFT, simple) {
    using namespace imajuscule;
    using namespace imajuscule::fft;
    
    constexpr auto N = 8;
    FFTVec<double> res, input{{
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
    auto roots = compute_roots_of_unity<double>(N);
    Algo<double> a(roots);
    a.run(input.begin(), res.begin(), N, 1 );
 
    // testing against ref found here : https://rosettacode.org/wiki/Fast_Fourier_transform
    
    ASSERT_NEAR(0, res[0].imag(), 1e-12);
    ASSERT_NEAR(0, res[2].imag(), 1e-12);
    ASSERT_NEAR(0, res[2].real(), 1e-12);
    ASSERT_NEAR(0, res[4].imag(), 1e-12);
    ASSERT_NEAR(0, res[4].real(), 1e-12);
    ASSERT_NEAR(0, res[6].imag(), 1e-12);
    ASSERT_NEAR(0, res[6].real(), 1e-12);
    
    ASSERT_NEAR(4, res[0].real(), 1e-12);
    ASSERT_NEAR(1, res[1].real(), 1e-12);
    ASSERT_NEAR(1, res[3].real(), 1e-12);
    ASSERT_NEAR(1, res[5].real(), 1e-12);
    ASSERT_NEAR(1, res[7].real(), 1e-12);
    ASSERT_NEAR(-2.41421, res[1].imag(), 1e-4);
    ASSERT_NEAR( 2.41421, res[7].imag(), 1e-4);
    ASSERT_NEAR(-0.414214, res[3].imag(), 1e-4);
    ASSERT_NEAR( 0.414214, res[5].imag(), 1e-4);
}

TEST(FFT, big) {
    using namespace imajuscule;
    using namespace imajuscule::fft;
    
    constexpr auto N = 4096;
    FFTVec<double> res, input;
    
    res.resize(N);
    input.resize(N);
    auto roots = compute_roots_of_unity<double>(N);
    Algo<double> a(roots);
    a.run(input.begin(), res.begin(), N, 1 );
}

