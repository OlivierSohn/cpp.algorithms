
namespace imajuscule {
    namespace testdspconv {
        
        template<typename Convolution>
        void testDirac() {
            using namespace fft;
            
            using T = typename Convolution::FPT;
            
            const std::vector<T> coefficients {{ .9,.8,.7,.6,.3,.2,.1 }};

            for(int i=0; i<5;i++){
                auto const part_size = pow2(i);
                Convolution conv;
                
                conv.set_partition_size(part_size);
                conv.setCoefficients(coefficients);
                
                // feed a dirac
                conv.step(1);
                for(int i=1; i<conv.getComputationPeriodicity(); ++i) {
                    conv.step(0);
                }
                std::vector<T> results;
                for(int i=0; i<coefficients.size(); ++i) {
                    results.push_back(conv.get());
                    conv.step(0);
                }
                
                // we sum the results of each partition so the epsilon is (worst case) :
                auto eps = conv.countPartitions() * getFFTEpsilon<T>(conv.get_fft_length());
                ASSERT_EQ(coefficients.size(), results.size());
                for(auto j=0; j<results.size(); ++j) {
                    ASSERT_NEAR(coefficients[j], results[j], eps) << j;
                }
            }
        }

        void testDirac() {
            //testDirac<FFTConvolution<float>>();
            //testDirac<FFTConvolution<double>>();
            testDirac<PartitionnedFFTConvolution<float>>();
            testDirac<PartitionnedFFTConvolution<double>>();
        }
}
}

TEST(Convolution, dirac) {
    using namespace imajuscule;
    using namespace imajuscule::testdspconv;
    
    testDirac();
}
