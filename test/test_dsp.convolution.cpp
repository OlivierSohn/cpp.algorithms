
namespace imajuscule::audio {
  namespace testdspconv {
    
    
    template<typename T>
    struct ConvolutionTraits {
      static constexpr bool supportsOddCountOfCoefficients = true;
      static constexpr bool evenIndexesAreApproximated = false;
      
      template <typename F>
      static void adaptCoefficients(F & v) {
      }
    };
    
    template<LatencySemantic L, typename T>
    struct ConvolutionTraits<SubSampled<L, T>> {
      static constexpr bool supportsOddCountOfCoefficients = false;
      static constexpr bool evenIndexesAreApproximated = true;
      
    };
    
#ifdef NDEBUG
    constexpr auto end_index = 16;
#else
      constexpr auto end_index = 15;
#endif

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
        case 14: return mkTestCoeffs<T>(200);
          // this should exceed the GPU capacity to do the fft in one kernel only,
          // and force partitionning:
        case 15: return mkTestCoeffs<T>(20000);
      }
      throw std::logic_error("coeff index too big");
    }
      
    template<typename Convolution, typename Coeffs, typename Input, typename Output>
    void testGeneric(Convolution & conv, Coeffs const & coefficients, Input const & input, Output const & expectedOutput, int vectorLength) {
      using T = typename Convolution::FPT;
      using namespace fft;

      if(!conv.isValid())
      {
        /*
         std::cout << std::endl << "Not testing invalid setup for "; COUT_TYPE(Convolution);
         std::cout << std::endl <<
         "coefficient size : " << coefficients.size() << std::endl <<
         "partition size : " << conv.getBlockSize() << std::endl <<
         "partition count : " << conv.countPartitions() << std::endl;
         */
        return;
      }
      
      using Tr = ConvolutionTraits<Convolution>;
      if(!Tr::supportsOddCountOfCoefficients
         && 1 == coefficients.size() % 2) {
        return;
      }
      conv.setCoefficients(coefficients);

      if(!conv.isValid())
      {
        // can happen for partitionned convolution on the cpu, when the count of coefficients is too low
        return;
      }
      

      
      std::vector<T> output;
      output.reserve(expectedOutput.size());

      auto const eps = conv.getEpsilon();
      int idxStep=0;
      auto step = [&idxStep, &input]() {
        return (idxStep < input.size()) ? input[idxStep] : ((T)0.);
      };
      
      for(; idxStep<conv.getLatency().toInteger(); ++idxStep) {
        auto res = conv.step(step());
        if(std::abs(res) > 1000*eps) {
          LG(INFO,"");
        }
        ASSERT_NEAR(0.f, res, 1000*eps); // assumes that no previous signal has been fed
      }
        std::vector<T> inputVec;
        inputVec.reserve(expectedOutput.size());
        for(;inputVec.size() != expectedOutput.size();++idxStep) {
            inputVec.push_back(step());
        }
        if(vectorLength) {
            // initialize with zeros
            output.resize(inputVec.size(), {});
            
            // then add
            int const nFrames = inputVec.size();
            for(int i=0; i<nFrames; i += vectorLength) {
                conv.stepAddVectorized(inputVec.data()+i,
                                      output.data()+i,
                                      std::min(vectorLength, nFrames-i));
            }
        }
        else {
            for(T i : inputVec) {
              output.push_back(conv.step(i));
            }
        }
      
      using Tr = ConvolutionTraits<Convolution>;
      for(auto j=0; j<output.size(); ++j) {
        if(Tr::evenIndexesAreApproximated && (0 == j%2)) {
          continue;
        }
        if(!areNear(expectedOutput[j], output[j], 1000*eps)) { // uses relative error
          std::cout << std::endl << "... "; COUT_TYPE(Convolution);
          std::cout << std::endl << "coefficient size : " << coefficients.size() << std::endl;
          ASSERT_NEAR(expectedOutput[j], output[j], 1000*eps); // doesn't use relative error, only absolute.
        }
      }
    }
    
    template<typename Convolution>
    void testDirac2(int coeffs_index, Convolution & conv) {
      
      if(!conv.isValid()) {
        return;
      }
      
      using T = typename Convolution::FPT;
      
      test(conv, makeCoefficients<T>(coeffs_index));
    }
    
    template<typename T, template<typename> typename Allocator, typename Tag>
    void testDirac() {
      using namespace fft;
      
        int end = end_index;
        if constexpr (!std::is_same_v<Tag, fft::Fastest>) {
            --end;
        }
      for(int i=0; i<end; ++i) {
          LG(INFO,"index %d", i);
        int const countCoeffs = makeCoefficients<T>(i).size();
          /*
          {
            auto c = FFTConvolution<T, Allocator, Tag>{};
            c.setup({static_cast<int>(ceil_power_of_two(countCoeffs))});
            testDirac2(i, c);
          }
           */
/*
          {
              for(int firstSz=1; firstSz<32; firstSz *= 2) {
                  auto c = CustomScaleConvolution<FFTConvolutionCore<T, Allocator, Tag>>{};
                  using ScalingParam = typename decltype(c)::SetupParam::ScalingParam;
                  std::vector<ScalingParam> params;
                  int remainingCoeffs = countCoeffs;
                  int sz = firstSz;
                  while(remainingCoeffs > 0) {
                      ScalingParam s {
                          // cannot do std::min(sz, remainingCoeffs) here else the latency would be wrong
                          sz,
                          {sz}
                      };
                      
                      params.push_back(s);
                      
                      remainingCoeffs -= sz;
                      sz *= 2;
                  }
                  c.setup({params});
                  testDirac2(i, c);
              }
          }
          // same as above but with some PartitionnedFFTConvolutionCRTP
          {
              for(int firstSz=1; firstSz<32; firstSz *= 2) {
                  auto c = CustomScaleConvolution<FFTConvolutionIntermediate < PartitionnedFFTConvolutionCRTP<T, Allocator, Tag> >>{};
                  using ScalingParam = typename decltype(c)::SetupParam::ScalingParam;
                  std::vector<ScalingParam> params;
                  int remainingCoeffs = countCoeffs;
                  int sz = firstSz;
                  while(remainingCoeffs > 0) {
                      int countCoeffs = std::min(sz, remainingCoeffs); // last coefficient group may be padded with 0s
                      ScalingParam s {
                          countCoeffs,
                          {
                              sz,
                              countPartitions(countCoeffs, sz)
                          }
                      };
                      
                      params.push_back(s);
                      
                      remainingCoeffs -= sz;
                      sz *= 2;
                  }
                  c.setup({params});
                  testDirac2(i, c);
              }
          }
          // same as above but skipping one scale out of 2
          {
              for(int firstSz=1; firstSz<32; firstSz *= 2) {
                  auto c = CustomScaleConvolution<FFTConvolutionIntermediate < PartitionnedFFTConvolutionCRTP<T, Allocator, Tag> >>{};
                  using ScalingParam = typename decltype(c)::SetupParam::ScalingParam;
                  std::vector<ScalingParam> params;
                  int remainingCoeffs = countCoeffs;
                  int sz = firstSz;
                  while(remainingCoeffs > 0) {
                      int const countCoeffs = std::min(sz*3, remainingCoeffs); // last coefficient group may be padded with 0s
                      ScalingParam s {
                          countCoeffs,
                          {
                              sz,
                              countPartitions(countCoeffs, sz)
                          }
                      };

                      params.push_back(s);

                      remainingCoeffs -= countCoeffs;
                      sz *= 4;
                  }
                  c.setup({params});
                  testDirac2(i, c);
              }
          }
        {
          auto c = AsyncCPUConvolution<FIRFilter<T, Allocator>, PolicyOnWorkerTooSlow::Wait>{};
            std::vector<int> submisionPeriods{
                -1, // invalid
                0, // invalid
                1,
                10,
                100
            };
          std::vector<int> queueSizes{
              -1, // invalid
              0, // invalid
              1,
              10,
              100
          };
            for(auto submisionPeriod:submisionPeriods) {
                for(auto queueSize:queueSizes) {
                    c.setup({
                        submisionPeriod,
                        queueSize,
                        {countCoeffs}
                    });
                    testDirac2(i, c);
                }
            }
        }*/
          /*
        if(i<10)
        {
          auto c = Delayed<FIRFilter<T, Allocator>>{};
          c.setup({10,{countCoeffs}});
          testDirac2(i, c);
        }
           */
        /*
        // TODO test T=double on a machine that has support for doubles on the gpu (cl_khr_fp64)
        {
          auto c = FIRFilterGPU<T>{};
          testDirac2(i, c);
        }
        {
          auto c = FIRFilterGPUAsync<T>{};
          testDirac2(i, c);
        }
        {
          auto c = FIRFilterGPUAsyncN<T>{};
          c.setup({10});
          testDirac2(i, c);
        }
         */
        // TODO test this on a machine that has support for doubles on the gpu (cl_khr_fp64)
        {
            // clBuildProgram is quite slow so I deactivate this
            /*
          auto c = PartitionnedFIRFilterGPUAsyncN<T>{};
          c.setup({10});
          testDirac2(i, c);
             */
        }
          /*
        {
            // min partition size of 4 else it is not valid (countGrains() < blockSize())
            int const szPartition = std::max(4,
                                             static_cast<int>(floor_power_of_two(countCoeffs/10)));
            int const latencyLateHandler = 2*szPartition - 1;
            int const nLateCoeffs = std::max(0, countCoeffs - latencyLateHandler);
            int const nEarlyCoeffs = countCoeffs - nLateCoeffs;
            auto c = mkRealTimeConvolutionSubsampled<T, Allocator, Tag>(countCoeffs,
                                                                        mkNaiveScaling(1, nEarlyCoeffs),
                                                                        szPartition,
                                                                        nLateCoeffs);
            testDirac2(i, c);
        }*/
      }
    }
  }
}

TEST(Convolution, dirac) {
  using namespace imajuscule;
  using namespace imajuscule::audio::testdspconv;
  
  for_each(fft::Tags, [](auto t) {
    testDirac<float, a64::Alloc, decltype(t)>();
    testDirac<double, a64::Alloc, decltype(t)>();
  });
}

TEST(Convolution, testComputeQueueSize) {
    using namespace imajuscule::audio;
    {
        int resultZeroIteration = computeQueueSize([](){ return 1.; }, 1., 0);
        ASSERT_EQ(1, resultZeroIteration);
    }
    {
        int resultSameSpeed = computeQueueSize([](){ return 12.; }, 12., 100);
        ASSERT_EQ(1, resultSameSpeed);
    }
    {
        int resultProcessingFaster = computeQueueSize([](){ return 1.; }, 12., 100);
        ASSERT_EQ(1, resultProcessingFaster);
    }
    {
        int i=-1;
        int resultProcessingFaster = computeQueueSize([&i]() mutable {
            ++i;
            switch(i%10) {
                case 3:
                    return 2.5;
                default:
                    return 0.01;
            }
            return 1.;
        }, 1., 100);
        ASSERT_EQ(3, resultProcessingFaster);
    }
    {
        int i=-1;
        int resultProcessingFaster = computeQueueSize([&i]() mutable {
            ++i;
            switch(i%10) {
                case 1:
                case 2:
                    return 0.01;
                case 3:
                    return 2.5;
                case 4:
                case 5:
                case 6:
                    return 0.01;
                case 7:
                    return 2.5;
                default:
                    return 0.01;
            }
            return 1.;
        }, 1., 100);
        ASSERT_EQ(3, resultProcessingFaster);
    }
    {
        int i=-1;
        int resultProcessingFaster = computeQueueSize([&i]() mutable {
            ++i;
            switch(i%10) {
                case 1:
                case 2:
                    return 0.01;
                case 3:
                    return 2.5;
                case 4:
                    return 3.0;
                case 5:
                case 6:
                case 7:
                    return 0.01;
                default:
                    return 0.01;
            }
            return 1.;
        }, 1., 100);
        ASSERT_EQ(5, resultProcessingFaster);
    }
    {
        int i=-1;
        int resultProcessingFaster = computeQueueSize([&i]() mutable {
            ++i;
            switch(i%10) {
                case 1:
                case 2:
                    return 0.01;
                case 3:
                    return 2.5;
                case 4:
                    return 0.01;
                case 5:
                    return 3.0;
                case 6:
                case 7:
                    return 0.01;
                default:
                    return 0.01;
            }
            return 1.;
        }, 1., 100);
        ASSERT_EQ(4, resultProcessingFaster);
    }
}
