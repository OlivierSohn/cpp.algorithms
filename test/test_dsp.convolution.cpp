
namespace imajuscule {
    namespace testdspconv {
        
        template<typename T>
        auto makeNonTrivialCoefficients() {
            return cacheline_aligned_allocated::vector<T>{{ .9,.8,.7,.6,.3,.2,.1,0. }};
        }
        
        template<typename T>
        auto makeTrivialCoefficients() {
            return cacheline_aligned_allocated::vector<T>{{ 1. }};
        }
        
        template<typename T>
        auto makeTrivialCoefficients2() {
            return cacheline_aligned_allocated::vector<T>{{ -1. }}; // makes partitionned test fail
        }
        
        template<typename T>
        auto makeCoefficients(int coeffs_index) {
            switch(coeffs_index) {
                case 0: return makeTrivialCoefficients<T>();
                case 1: return makeTrivialCoefficients2<T>();
                default:
                    return makeNonTrivialCoefficients<T>();
            }
        }
        
        template<typename Convolution, typename Coeffs>
        bool test(Convolution & conv, Coeffs const & coefficients) {
            using T = typename Convolution::FPT;
            using namespace fft;
            
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
            EXPECT_EQ(coefficients.size(), results.size());
            for(auto j=0; j<results.size(); ++j) {
                EXPECT_NEAR(coefficients[j], results[j], eps);
            }
            return false;
        }
        
        template<typename Convolution>
        bool testDiracPartitionned(int coeffs_index) {
            
            using T = typename Convolution::FPT;
            
            const auto coefficients = makeCoefficients<T>(coeffs_index);
            
            int i=0;
            for(int i=0; i<5;i++)
            {
                auto const part_size = pow2(i);
                Convolution conv;
                
                conv.set_partition_size(part_size);
                test(conv, coefficients);
            }
            return false;
        }
        
        template<typename Convolution>
        bool testDirac2(int coeffs_index) {
            
            using T = typename Convolution::FPT;
            
            const auto coefficients = makeCoefficients<T>(coeffs_index);
            Convolution conv;
            
            test(conv, coefficients);
            return false;
        }
        
        template<typename Tag>
        bool testDirac() {
            using namespace fft;
            for(int i=0; i<3; ++i) {
                testDirac2<FFTConvolution<float, Tag>>(i);
                testDirac2<FFTConvolution<double, Tag>>(i);
                testDiracPartitionned<PartitionnedFFTConvolution<float, Tag>>(i);
                testDiracPartitionned<PartitionnedFFTConvolution<double, Tag>>(i);
            }
            return false;
        }
    }
}

TEST(Convolution, dirac) {
    using namespace imajuscule;
    using namespace imajuscule::testdspconv;
    
    for_each(fft::Tags, [](auto t) {
        testDirac<decltype(t)>();
    });
}

