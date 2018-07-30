
namespace imajuscule {
    namespace testspatialize {
        constexpr auto end_index = 4;
        
        template<typename T>
        a64::vector<T> makeCoefficients(int coeffs_index) {
            switch(coeffs_index) {
                case 0: return {{ +1. }};
                case 1: return {{ -1. }};
                case 2: return {{ .9,.8,.7,.6,.3,.2,.1,0. }};
                case 3: {
                    constexpr auto sz = 2000;
                    a64::vector<T> v(sz);
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
        
        template<typename Spatialize, typename T = typename Spatialize::FPT>
        void test(Spatialize & spatialize, std::array<a64::vector<T>, Spatialize::nEars> const & ear_signals) {
            using namespace fft;
            constexpr auto nEars = Spatialize::nEars;
            
            spatialize.setMaxMultiplicationGroupLength();

            if(!spatialize.isValid()) {
                /*
                 std::cout << std::endl << "Not testing invalid setup for "; COUT_TYPE(Convolution);
                std::cout << std::endl << 
                "coefficient size : " << coefficients.size() << std::endl <<
                "partition size : " << conv.getBlockSize() << std::endl <<
                "partition count : " << conv.countPartitions() << std::endl;
                 */
                return;
            }
            
            // feed a dirac
            std::vector<T> input;
            input.resize(spatialize.countSources());
            std::fill(input.begin(), input.end(), 1);
            
            spatialize.step(input.data());
            for(int i=0; i<spatialize.getLatency(); ++i) {
                std::fill(input.begin(), input.end(), 0);
                spatialize.step(input.data());
            }
            
            std::fill(input.begin(), input.end(), 0);

            std::vector<std::array<T, nEars>> results;
            for(int i=0; i<ear_signals[0].size(); ++i) {
                results.emplace_back();
                spatialize.get(results.back().data());
                spatialize.step(input.data());
            }

            auto eps = spatialize.getEpsilon();
            ASSERT_EQ(ear_signals[0].size(), results.size());
            for(auto j=0; j<results.size(); ++j) {
                int i=0;
                for(auto r : results[j]) {
                    ASSERT_NEAR(ear_signals[i][j], r, eps);
                    ++i;
                }
            }
        }
        
        template<typename T, typename F>
        void testSpatialized(int coeffs_index, F f) {
            const auto coefficients = makeCoefficients<T>(coeffs_index);
            
            if(coefficients.size() < 1024) {
                for(int i=0; i<5;i++)
                {
                    const auto part_size = pow2(i);
                    f(part_size, coefficients);
                }
            }
            else {
                constexpr auto part_size = 256;
                f(part_size, coefficients);
            }
        }
        
        template<typename Convolution>
        void testDiracFinegrainedPartitionned(int coeffs_index) {
            
            auto f = [](int part_size, auto const & coefficients)
            {
                typename std::remove_const<typename std::remove_reference<decltype(coefficients)>::type>::type zero_coeffs, opposite_coeffs;
                zero_coeffs.resize(coefficients.size());

                {
                    constexpr auto nOutMono = 1;
                    audio::Spatializer<nOutMono, Convolution> spatialized;
                    
                    spatialized.set_partition_size(part_size);
                    
                    spatialized.addSourceLocation({{coefficients}});
                    
                    test(spatialized, {{coefficients}});
                }
                {
                    constexpr auto nOutStereo = 2;
                    audio::Spatializer<nOutStereo, Convolution> spatialized;
                    
                    spatialized.set_partition_size(part_size);
                    
                    // no cross-talk :
                    // source1 has only left component
                    // source2 has only right component
                    spatialized.addSourceLocation({{coefficients, zero_coeffs}});
                    spatialized.addSourceLocation({{zero_coeffs, coefficients}});
                    
                    test(spatialized, {{coefficients, coefficients}});
                }
                
                {
                    constexpr auto nOutStereo = 2;
                    audio::Spatializer<nOutStereo, Convolution> spatialized;
                    
                    spatialized.set_partition_size(part_size);

                    opposite_coeffs = coefficients;
                    for(auto & o : opposite_coeffs) {
                        o = -o;
                    }
                    // extreme cross-talk :
                    // source1 right component is the opposite of source 2
                    // source2 right component is the opposite of source 1
                    spatialized.addSourceLocation({{coefficients, opposite_coeffs}});
                    spatialized.addSourceLocation({{opposite_coeffs, coefficients}});
                    
                    test(spatialized, {{zero_coeffs,zero_coeffs}});
                }
                {
                    constexpr auto nOutStereo = 2;
                    audio::Spatializer<nOutStereo, Convolution> spatialized;
                    
                    spatialized.set_partition_size(part_size);
                    opposite_coeffs = coefficients;
                    for(auto & o : opposite_coeffs) {
                        o = -o;
                    }

                    spatialized.addSourceLocation({{coefficients, opposite_coeffs}});
                    spatialized.addSourceLocation({{zero_coeffs, coefficients}});
                    
                    test(spatialized, {{coefficients,zero_coeffs}});
                }
            };
            
            testSpatialized<typename Convolution::FPT>(coeffs_index, f);
        }
        
        template<typename Tag>
        void testDirac() {
            using namespace fft;
            for(int i=0; i<end_index; ++i) {
                testDiracFinegrainedPartitionned<FinegrainedPartitionnedFFTConvolution<float, Tag>>(i);
                testDiracFinegrainedPartitionned<FinegrainedPartitionnedFFTConvolution<double, Tag>>(i);
            }
        }
    }
}

TEST(Spatialization, dirac) {
    using namespace imajuscule;
    using namespace imajuscule::testspatialize;
    
    for_each(fft::Tags, [](auto t) {
        testDirac<decltype(t)>();
    });
}

