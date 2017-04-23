
namespace imajuscule {
    namespace testdspconv {
        constexpr auto end_index = 4;
        
        template<typename T>
        a_64::vector<T> makeCoefficients(int coeffs_index) {
            switch(coeffs_index) {
                case 0: return {{ +1. }};
                case 1: return {{ -1. }};
                case 2: return {{ .9,.8,.7,.6,.3,.2,.1,0. }};
                case 3: {
                    constexpr auto sz = 80000;
                    a_64::vector<T> v(sz);
                    auto index = 0;
                    for(auto & value: v) {
                        value = (sz - index) / static_cast<T>(sz);
                        ++index;
                    }
                    return std::move(v);
                }
            }
            throw std::logic_error("coeff index too big");
        }
        
        template<typename Convolution, typename Coeffs>
        bool test(Convolution & conv, Coeffs const & coefficients) {
            using T = typename Convolution::FPT;
            using namespace fft;
            
            conv.setCoefficients(coefficients);
            
            // feed a dirac
            conv.step(1);
            for(int i=1; i<conv.getLatency(); ++i) {
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
            
            if(coefficients.size() < 1024) {
                for(int i=0; i<5;i++)
                {
                    auto const part_size = pow2(i);
                    Convolution conv;
                    
                    conv.set_partition_size(part_size);
                    test(conv, coefficients);
                }
            }
            else {
                Convolution conv;
                
                conv.set_partition_size(1024);
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
            for(int i=0; i<end_index; ++i) {
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

