
namespace imajuscule {
    namespace testdspconv {
        constexpr auto end_index = 15;

        template<typename T>
        a64::vector<T> makeCoefficients(int coeffs_index) {
            switch(coeffs_index) {
                case 0: return {{ +1. }};
                case 1: return {{ -1. }};
              case 2: return {{ .9, }};
              case 3: return {{ .9,.8 }};
              case 4: return {{ .9,.8,.7 }};
              case 5: return {{ .9,.8,.7,.6 }};
              case 6: return {{ .9,.8,.7,.6,.5 }};
              case 7: return {{ .9,.8,.7,.6,.5,.4 }};
              case 8: return {{ .9,.8,.7,.6,.5,.4,.3 }};
              case 9: return {{ .9,.8,.7,.6,.5,.4,.3,.2 }};
              case 10: return {{ .9,.8,.7,.6,.5,.4,.3,.2,.1 }};
              case 11: return {{ .9,.8,.7,.6,.5,.4,.3,.2,.1,.05 }};
              case 12: return {{ .9,.8,.7,.6,.5,.4,.3,.2,.1,.05,.025 }};
              case 13: return {{ .9,.8,.7,.6,.5,.4,.3,.2,.1,.05,.025,.01 }};
                case 14: {
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

        template<typename Convolution, typename Coeffs>
        void test(Convolution & conv, Coeffs const & coefficients) {
            using T = typename Convolution::FPT;
            using namespace fft;

            if(!conv.isValid()) {
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
            conv.step(1);
            for(int i=0; i<conv.getLatency(); ++i) {
                conv.step(0);
            }

            std::vector<T> results;
            for(int i=0; i<coefficients.size(); ++i) {
                results.push_back(conv.get());
                conv.step(0);
            }

            auto eps = conv.getEpsilon();
            ASSERT_EQ(coefficients.size(), results.size());
            for(auto j=0; j<results.size(); ++j) {
                ASSERT_NEAR(coefficients[j], results[j], eps);
            }
        }

        template<typename T, typename F>
        void testPartitionned(int coeffs_index, F f) {
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

        enum class TestFinegrained {
            Begin,

            Low = Begin,
            Med,
            High,

            End
        };

        template<typename Convolution>
        void testDiracFinegrainedPartitionned(int coeffs_index) {

            auto f = [](int part_size, auto const & coefficients)
            {
                for(auto type = TestFinegrained::Begin;
                    type != TestFinegrained::End;
                    increment(type))
                {
                    Convolution conv;

                    conv.set_partition_size(part_size);
                    conv.setCoefficients(coefficients);

                    range<int> r {
                        conv.getLowestValidMultiplicationsGroupSize(),
                        conv.getHighestValidMultiplicationsGroupSize()
                    };

                    switch(type) {
                        case TestFinegrained::Low:
                            conv.setMultiplicationGroupLength(r.getMin());
                            break;
                        case TestFinegrained::High:
                            conv.setMultiplicationGroupLength(r.getMax());
                            break;
                        case TestFinegrained::Med:
                            conv.setMultiplicationGroupLength(r.getExpCenter());
                            break;
                        default:
                            throw std::logic_error("not supported");
                    }
                    test(conv, coefficients);
                }
            };

            testPartitionned<typename Convolution::FPT>(coeffs_index, f);
        }

        template<typename Convolution>
        void testDiracPartitionned(int coeffs_index) {

            auto f = [](int part_size, auto const & coefficients){
                Convolution conv;

                conv.set_partition_size(part_size);
                conv.setCoefficients(coefficients);
                test(conv, coefficients);
            };

            testPartitionned<typename Convolution::FPT>(coeffs_index, f);
        }

        template<typename Convolution>
        void testDirac2(int coeffs_index, Convolution & conv) {

            using T = typename Convolution::FPT;

            const auto coefficients = makeCoefficients<T>(coeffs_index);
            conv.setCoefficients(coefficients);

            test(conv, coefficients);
        }
      
        template<typename T, typename FFTTag = fft::Fastest>
        auto mkRealTimeConvolution() {
          auto c = ZeroLatencyBruteFineGrainedPartitionnedConvolution<T, FFTTag>{};
          c.set_partition_size(4);
          c.applySetup(FinegrainedSetupParam{1,0});
          return c;
        }

      
      template<typename T, typename FFTTag = fft::Fastest>
      auto mkRealTimeConvolution2() {
        auto c = ZeroLatencyScaledFineGrainedPartitionnedConvolution<T, FFTTag>{};
        c.set_partition_size(4);
        c.applySetup(FinegrainedSetupParam{1,0});
        return c;
      }

        template<typename Tag>
        bool testDirac() {
            using namespace fft;
            for(int i=0; i<end_index; ++i) {
                testDiracFinegrainedPartitionned<FinegrainedPartitionnedFFTConvolution<float, Tag>>(i);
                testDiracFinegrainedPartitionned<FinegrainedPartitionnedFFTConvolution<double, Tag>>(i);
              {
                auto c = ZeroLatencyBruteConvolution<float>{};
                testDirac2(i, c);
              }
              {
                auto c = ZeroLatencyBruteConvolution<double>{};
                testDirac2(i, c);
              }
              {
                auto c = FFTConvolution<float, Tag>{};
                testDirac2(i, c);
              }
              {
                auto c = FFTConvolution<double, Tag>{};
                testDirac2(i, c);
              }
              {
                auto c = mkRealTimeConvolution<float, Tag>();
                testDirac2(i, c);
              }
              {
                auto c = mkRealTimeConvolution<double, Tag>();
                testDirac2(i, c);
              }
              LG(INFO,"2 ?");
              {
                auto c = mkRealTimeConvolution2<float, Tag>();
                testDirac2(i, c);
              }
              {
                auto c = mkRealTimeConvolution2<double, Tag>();
                testDirac2(i, c);
              }
              LG(INFO,"2");
                testDiracPartitionned<PartitionnedFFTConvolution<float, Tag>>(i);
                testDiracPartitionned<PartitionnedFFTConvolution<double, Tag>>(i);
                //testDiracPartitionned<ScalingPartitionnedFFTConvolution<float, Tag>>(i);
                //testDiracPartitionned<ScalingPartitionnedFFTConvolution<double, Tag>>(i);
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
