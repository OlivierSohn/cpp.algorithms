

namespace imajuscule::audio {
static void testByCbSize(int dirac_offset, int cb_change) {
    constexpr int nAudioOut = 1;
    
    // on choisit un type de reverb pour lequel l'optimisation est rapide
    ConvReverbsByBlockSize<Reverbs<nAudioOut, ReverbType::Offline, PolicyOnWorkerTooSlow::Wait>> r;

    ASSERT_EQ(r.getWetRatioWithoutStepping(), 1.);
    r.abruptlySetConvolutionReverbWetRatio(0.5);
    ASSERT_EQ(r.getWetRatioWithoutStepping(), 0.5);
    r.abruptlySetConvolutionReverbWetRatio(1.);
    ASSERT_EQ(r.getWetRatioWithoutStepping(), 1.);

    
    const int nSources = 1;
    auto coeffs = mkCoefficientsTriangle(10000);
    std::map<int, ConvReverbOptimizationReport> m;
    
    ASSERT_TRUE(coeffs[0] == 1.);
    
    XFFTsCostsFactors factors;
    r.setConvolutionReverbIR(nSources,
                           DeinterlacedBuffers<double>(coeffs, 1),
                           1024,
                           0,
                           44100.,
                           m,
                           factors
                           );
    ASSERT_EQ(0, m.count(2048));
    ASSERT_EQ(1, m.count(1024));
    ASSERT_EQ(1, m.count(512));
    ASSERT_EQ(1, m.count(256));
    ASSERT_EQ(1, m.count(4));
    ASSERT_EQ(1, m.count(2));
    ASSERT_EQ(1, m.count(1));
    
    ASSERT_EQ(r.getWetRatioWithoutStepping(), 1.);
    r.abruptlySetConvolutionReverbWetRatio(0.5);
    ASSERT_EQ(r.getWetRatioWithoutStepping(), 0.5);
    r.abruptlySetConvolutionReverbWetRatio(1.);
    ASSERT_EQ(r.getWetRatioWithoutStepping(), 1.);

    std::vector<double> io_buffer;
    io_buffer.push_back(1.);
    std::vector<double*> io_buffers;
    io_buffers.push_back(io_buffer.data());

    std::vector<double> work_buffer;
    work_buffer.resize(1);
    std::vector<double*> work_buffers;
    work_buffers.push_back(work_buffer.data());

    ASSERT_EQ(1., io_buffer[0]);

    // no block size has been declared, hence only zeros should be output:
    r.apply(io_buffers.data(),
            1,
            work_buffers.data(),
            1,
            1);
    
    ASSERT_EQ(0., io_buffer[0]);
    
    int sz = r.declareBlockSize(1024);
    
    ASSERT_EQ(1024, sz);
    
    std::vector<double> input; // offset dirac
    std::vector<double> expected_output;

    for(int i=0; i<dirac_offset; ++i) {
        input.push_back(0.);
        expected_output.push_back(0.);
    }
    input.push_back(1.);
    expected_output.push_back(coeffs[0]);
    for(int i=1; i<2*coeffs.size(); ++i) {
        input.push_back(0.);
        expected_output.push_back(i<coeffs.size() ? coeffs[i] : 0.);
    }

    for(int i=0; i<input.size(); ++i) {
        
        if(i==cb_change) {
            sz = r.declareBlockSize(128);
            ASSERT_EQ(128, sz);
        }
        
        io_buffer[0] = input[i];
        
        r.apply(io_buffers.data(),
                1,
                work_buffers.data(),
                1,
                1);
        
        if(!areNear(expected_output[i], io_buffer[0], 1e-9)) {
            ASSERT_NEAR(expected_output[i], io_buffer[0], 1e-9);
        }
    }

}
} // NS imajuscule

TEST(ReverbsByCbSize, dirac) {
    using namespace imajuscule::audio;
    
    for(int dirac_offset = 0; dirac_offset < 5; ++dirac_offset) {
        for(int cb_change = 0; cb_change < 5; ++cb_change) {
            testByCbSize(dirac_offset, cb_change);
        }
    }
}
