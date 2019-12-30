

static inline float randf(float high = 1.f, float low = 0.f)
{
  return low + ((high-low) * ((float)std::rand() / (float)(RAND_MAX + 1.f)));
}

namespace imajuscule {
  namespace testdspconv {
    
    template<typename It>
    void averageNeighbours(It it, It end) {
      while(it+1 < end) {
        *it = *(it+1) = 0.5f * (*it + *(it+1));
        it += 2;
      }
    }
    
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
        case 14: return mkTestCoeffs<T>(2000);
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
      
      for(; idxStep<conv.getLatency(); ++idxStep) {
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
        
  template<typename Convolution, typename Coeffs, typename Input, typename Output>
  void testGeneric(Convolution & conv, Coeffs const & coefficients, Input const & input, Output const & expectedOutput) {
      int sz = expectedOutput.size();
      
      std::set<int> vectorLengths;
      
      vectorLengths.insert(0);// no vectorization
      vectorLengths.insert(1);// minimal vectorization
      vectorLengths.insert(2);
      vectorLengths.insert(3);
      vectorLengths.insert(sz/5);
      vectorLengths.insert(sz/4);
      vectorLengths.insert(sz/4-1);
      vectorLengths.insert(sz/4+1);
      vectorLengths.insert(sz/2);
      vectorLengths.insert(sz);
      vectorLengths.insert(2*sz);

      
      for(auto vectorLength: vectorLengths) {
          if(vectorLength < 0) {
              continue;
          }
          testGeneric(conv, coefficients, input, expectedOutput, vectorLength);
      }
  }

    // we take 'coefficients' by value because we modify it in the function
    template<typename Convolution, typename Coeffs>
    void test(Convolution & conv, Coeffs const coefficients)
    {
      using T = typename Convolution::FPT;
      using Tr = ConvolutionTraits<Convolution>;

      // test using a dirac
      {
        auto diracInput = mkDirac<T>(coefficients.size());
        auto output = coefficients;
        Tr::adaptCoefficients(output);
        testGeneric(conv, coefficients, diracInput, output);
      }
      
      // test using a random (but reproducible) signal
      if(!Tr::evenIndexesAreApproximated)
      {
        // brute-force convolution is costly (N^2) for high number of coefficients, hence we use caching
        static std::map<Coeffs, std::pair<std::vector<T>, std::vector<T>>> cache;
        
        std::vector<T> randomInput, output;

        auto it = cache.find(coefficients);
        if(it == cache.end()) {
          const int sz = 5*coefficients.size();
          output.reserve(sz);
          randomInput.reserve(sz);
          std::srand(0); // to have reproducible random numbers
          for(int i=0; i<sz; ++i) {
            randomInput.push_back(randf(1.f, -1.f));
          }
          // produce expected output using brute force convolution
          {
            FIRFilter<T> filter;
            filter.setCoefficients(coefficients);
            EXPECT_EQ(0, filter.getLatency());
            for(int i=0; i<randomInput.size(); ++i) {
              output.push_back(filter.step(randomInput[i]));
            }
            // "flush" the reverb
            for(int i=0; i<coefficients.size(); ++i) {
              output.push_back(filter.step(0.));
            }
          }
          
          cache.emplace(coefficients, std::make_pair(randomInput, output));
        }
        else {
          randomInput = it->second.first;
          output = it->second.second;
        }
        
        Tr::adaptCoefficients(output);
        testGeneric(conv, coefficients, randomInput, output);
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
    
    template<typename T, typename Tag>
    void testDiracFinegrainedPartitionned(int coeffs_index) {
      
      auto f = [](int part_size, auto & coefficients)
      {
        for(auto type = TestFinegrained::Begin;
            type != TestFinegrained::End;
            increment(type))
        {
          FinegrainedPartitionnedFFTConvolution<T, Tag> conv;
          
          conv.setup({part_size, 1000, 0});
          if(!conv.isValid()) {
            continue;
          }
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
      
      testPartitionned<T>(coeffs_index, f);
    }
    
    template<typename T, typename Tag>
    void testDiracPartitionned(int coeffs_index) {
      
      auto f = [](int part_size, auto & coefficients){
        PartitionnedFFTConvolution<T,Tag> conv;
        
        conv.setup({part_size});

        test(conv, coefficients);
      };
      
      testPartitionned<T>(coeffs_index, f);
    }
    
    
    template<typename Convolution>
    void testDirac2(int coeffs_index, Convolution & conv) {
      
      if(!conv.isValid()) {
        return;
      }
      
      using T = typename Convolution::FPT;
      
      test(conv, makeCoefficients<T>(coeffs_index));
    }
    
    
    template<typename T, typename FFTTag = fft::Fastest>
    auto mkRealTimeConvolution(std::vector<Scaling> const & v, int partitionSize) {
      using C = ZeroLatencyScaledFineGrainedPartitionnedSubsampledConvolution<T, FFTTag>;
      using ScalingParam = typename C::SetupParam::AParam::BParam::ScalingParam;
      
      auto scalingParams = scalingsToParams<ScalingParam>(v);
      auto c = C{};
      c.setup(typename C::SetupParam
      {
        {
            {},
            {scalingParams}
        },
        {
          FinegrainedSetupParam{
              partitionSize, // partition size
              partitionSize*1000, // multiplication group size
              0 // phase
          },
          {
            0,
            {FinegrainedSetupParam{0,1,0},
              {
                0,
                {FinegrainedSetupParam{0,1,0},
                  {
                    0,
                    FinegrainedSetupParam{0,1,0}
                  }}
              }}
          }
        }
      }
                 );
      return c;
    }
    
    template<typename T, typename Tag>
    void testDirac() {
      using namespace fft;
      
      for(int i=0; i<end_index; ++i) {
          LG(INFO,"index %d", i);
        int const countCoeffs = makeCoefficients<T>(i).size();
        testDiracFinegrainedPartitionned<T, Tag>(i);
        if(i<10)
        {
          auto c = FIRFilter<T>{};
          testDirac2(i, c);
        }
          for(int dropped=0; dropped<5; ++dropped) {
              auto c = ScaleConvolution<FFTConvolutionCore<T, Tag>>{};
              c.setup({CountDroppedScales(dropped)});
              testDirac2(i, c);
          }
          // This is the exact same as above, but using CustomScaleConvolution
          {
              for(int firstSz=1; firstSz<32; firstSz *= 2) {
                  auto c = CustomScaleConvolution<FFTConvolutionCore<T, Tag>>{};
                  using ScalingParam = typename decltype(c)::SetupParam::ScalingParam;
                  std::vector<ScalingParam> params;
                  int remainingCoeffs = countCoeffs;
                  int sz = firstSz;
                  while(remainingCoeffs > 0) {
                      ScalingParam s;
                      // cannot do std::min(sz, remainingCoeffs) here else the latency would be wrong
                      s.countCoeffs = sz;
                      s.submissionPeriod = sz;
                      s.setupParam = {};
                      
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
                  auto c = CustomScaleConvolution<FFTConvolutionIntermediate < PartitionnedFFTConvolutionCRTP<T, Tag> >>{};
                  using ScalingParam = typename decltype(c)::SetupParam::ScalingParam;
                  std::vector<ScalingParam> params;
                  int remainingCoeffs = countCoeffs;
                  int sz = firstSz;
                  while(remainingCoeffs > 0) {
                      int countCoeffs = std::min(sz, remainingCoeffs); // last coefficient group may be padded with 0s
                      int submissionPeriod = sz;
                      ScalingParam s{countCoeffs,
                                     submissionPeriod,
                                     {submissionPeriod}};
                      
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
                  auto c = CustomScaleConvolution<FFTConvolutionIntermediate < PartitionnedFFTConvolutionCRTP<T, Tag> >>{};
                  using ScalingParam = typename decltype(c)::SetupParam::ScalingParam;
                  std::vector<ScalingParam> params;
                  int remainingCoeffs = countCoeffs;
                  int sz = firstSz;
                  while(remainingCoeffs > 0) {
                      int const countCoeffs = std::min(sz*3, remainingCoeffs); // last coefficient group may be padded with 0s
                      int const submissionPeriod = sz;
                      ScalingParam s{countCoeffs,
                                     submissionPeriod,
                                     {submissionPeriod}};
                      
                      params.push_back(s);

                      remainingCoeffs -= countCoeffs;
                      sz *= 4;
                  }
                  c.setup({params});
                  testDirac2(i, c);
              }
          }
        {
          auto c = AsyncCPUConvolution<FIRFilter<T>>{};
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
                        {}
                    });
                    testDirac2(i, c);
                }
            }
        }
        if(i<10)
        {
          auto c = Delayed<FIRFilter<T>>{};
          c.setup({10,{}});
          testDirac2(i, c);
        }
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
          auto c = PartitionnedFIRFilterGPUAsyncN<T>{};
          c.setup({10});
          testDirac2(i, c);
        }
        {
          auto c = FFTConvolution<T, Tag>{};
          testDirac2(i, c);
        }
        {
            // min partition size of 4 else it is not valid (countGrains() < blockSize())
            int const szPartition = std::max(4,
                                             static_cast<int>(floor_power_of_two(countCoeffs/10)));
            int const latencyLateHandler = 2*szPartition - 1;
            int const nLateCoeffs = std::max(0, countCoeffs - latencyLateHandler);
            int const nEarlyCoeffs = countCoeffs - nLateCoeffs;
            auto c = mkRealTimeConvolution<T, Tag>(mkNaiveScaling(1, nEarlyCoeffs), szPartition);
            testDirac2(i, c);
        }
        testDiracPartitionned<T, Tag>(i);
      }
    }
  }
}

TEST(Convolution, dirac) {
  using namespace imajuscule;
  using namespace imajuscule::testdspconv;
  
  for_each(fft::Tags, [](auto t) {
    testDirac<float, decltype(t)>();
    testDirac<double, decltype(t)>();
  });
}

// shows that the fft of the coefficients ** do not ** give the long term amplitude
// by frequency.
TEST(Convolution, freq) {
  using namespace imajuscule;
  using namespace imajuscule::fft;
  using namespace imajuscule::fft::slow_debug;
  using namespace imajuscule::testdspconv;
  
  //using Tag = Fastest;
  using Tag = imj::Tag;
  
  using ScopedContext = ScopedContext_<Tag, double>;
  using Algo = Algo_<Tag, double>;
  using RealSignal = typename fft::RealSignal_<Tag, double>::type;
  using CplxFreqs = typename fft::RealFBins_<Tag, double>::type;
  
  std::vector<Scaling> scalingParams{}; // leave it empty, we have only 8 coefficients
  auto c = mkRealTimeConvolution<double>(scalingParams, 4);
  constexpr auto N = 8;
  
  //a64::vector<double> coefficients{1., 0.707106, 0., -0.707106, -1., -0.707106, 0., 0.707106};
  //a64::vector<double> coefficients{1., 0.5, 0., -0.5, -1., -0.5, 0., 0.5};
  //a64::vector<double> coefficients{1., 0.75, 0.25, 0., -0.25, -0.5, -0.75, -1.0};
  a64::vector<double> coefficients{1., -0.5, 0.25, -0.125, 0.06, -0.03, 0.01, -0.005};
  /*
   corresponding norms:
   
   1.98
   1.35721
   0.894427
   0.714822
   0.66
   0.714822
   0.894427
   1.35721
   
   d     : 1.98    // constant
   f2 16 : 1.758
   f2 8  : 1.304
   f     : 1.304
   f2 6  : 0.9975
   h     : 0.8
   g     : 0.66
   */
  
  c.setCoefficients(coefficients);
  auto d = [](int i) {
    return 1.;
  };
  
  auto g = [](int i) {
    switch(i%2) {
      case 0: return 1.;
      case 1: return -1.;
    }
    Assert(0);
    return 0.;
  };
  
  auto h = [](int i) {
    switch(i%4) {
      case 0: return 1.;
      case 1: return 0.;
      case 2: return -1.;
      case 3: return 0.;
    }
    Assert(0);
    return 0.;
  };
  
  auto f = [](int i) {
    switch(i%8) {
      case 0: return 1.;
      case 1: return 0.707106;
      case 2: return 0.;
      case 3: return -0.707106;
      case 4: return -1.;
      case 5: return -0.707106;
      case 6: return 0.;
      case 7: return 0.707106;
    }
    Assert(0);
    return 0.;
  };
  
  auto f2 = [](double i, int j) {
    return cos(static_cast<double>(j) * 2 * M_PI * i / N);
  };
  
  auto f3 = [](int i, int j) {
    return (i%(2*j) >= j) ? 1.f : -1.f;
  };
  
  for(int j=1; j<20; ++j) {
    double maxOut = 0;
    for(int i=0; i<10000; ++i) {
      auto res = c.step(f3(i,j));
      if(i>1000) {
        if(maxOut < res)
          maxOut = res;
      }
    }
    LG(INFO,"%d : %f",j,maxOut);
  }
  
  ScopedContext setup(N);
  
  Algo fft_algo(setup.get());
  CplxFreqs fft_of_coeffs;
  fft_of_coeffs.resize(N);
  auto coeffVec = fft::RealSignal_<Tag, double>::make(coefficients);
  ASSERT_EQ(N, coeffVec.size());
  fft_algo.forward(coeffVec.begin(), fft_of_coeffs, N);
  auto unwrapped_fft_of_coeffs = unwrap_frequencies<Tag>(fft_of_coeffs, N);
  for(auto &e : unwrapped_fft_of_coeffs) {
    e *= 1 / Algo::scale;
  }
  for(auto const&e : unwrapped_fft_of_coeffs) {
    std::cout << abs(e) << std::endl;
  }
  LG(INFO,"%f", 0.f);
  
}


TEST(Convolution, testComputeQueueSize) {
    using namespace imajuscule;
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
