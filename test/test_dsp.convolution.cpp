template<typename T>
static inline std::vector<T> mkDirac(int sz, T amplitude = 1.) {
  std::vector<T> res;
  res.resize(sz);
  if(res.empty()) {
    throw std::logic_error("0-size");
  }
  res[0] = amplitude;
  return res;
}

static inline float randf(float high = 1.f, float low = 0.f)
{
  return low + ((high-low) * ((float)std::rand() / (float)(RAND_MAX + 1.f)));
}
template<typename T>
bool areNear(T a, T b, double eps) {
  if(std::abs(a) < eps && std::abs(b) < eps) {
    return true;
  }
  if(std::abs(a-b) < eps) {
    return true;
  }
  return std::abs(a-b)/std::abs(std::max(a,b)) < eps;
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
      
      template <typename F>
      static void adaptCoefficients(F & v) {
        averageNeighbours(v.begin(),v.end());
      }
    };
    
#ifdef NDEBUG
    constexpr auto end_index = 16;
#else
      constexpr auto end_index = 15;
#endif

    template<typename T>
    inline auto mkTestCoeffs(int const sz){
      a64::vector<T> v(sz);
      auto index = 0;
      for(auto & value: v) {
        value = (sz - index) / static_cast<T>(sz);
        ++index;
      }
      return std::move(v);
    }
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
    void testGeneric(Convolution & conv, Coeffs const & coefficients, Input const & input, Output const & expectedOutput)
    {
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
      int i=0;
      auto step = [&i, &input]() {
        return (i < input.size()) ? input[i] : ((T)0.);
      };
      
      for(; i<conv.getLatency(); ++i) {
        auto res = conv.step(step());
        if(std::abs(res) > 1000*eps) {
          LG(INFO,"");
        }
        ASSERT_NEAR(0.f, res, 1000*eps); // assumes that no previous signal has been fed
      }
      for(;output.size() != expectedOutput.size();++i) {
        output.push_back(conv.step(step()));
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
          
          applySetup(conv, {part_size, 1000, 0});
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
        
        applySetup(conv, {1.f,{part_size}});

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
    auto mkRealTimeConvolution() {
      using C = ZeroLatencyScaledFineGrainedPartitionnedConvolution<T, FFTTag>;
      auto c = C{};
      applySetup(c,typename C::SetupParam
      {
        {{},{}},
        {
          FinegrainedSetupParam{4,1000,0},
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
    
    template<typename Tag>
    bool testDirac() {
      using namespace fft;
      
      for(int i=0; i<end_index; ++i) {
          LG(INFO,"index %d", i);
        testDiracFinegrainedPartitionned<float, Tag>(i);
        testDiracFinegrainedPartitionned<double, Tag>(i);
        {
          auto c = FIRFilter<float>{};
          testDirac2(i, c);
        }
        {
          auto c = FIRFilter<double>{};
          testDirac2(i, c);
        }
        {
          auto c = Delayed<FIRFilter<double>>{};
          applySetup(c,{10,{}});
          testDirac2(i, c);
        }
        /*
        {
          auto c = FIRFilterGPU<float>{};
          testDirac2(i, c);
        }
        {
          // TODO test this on a machine that has support for doubles on the gpu (cl_khr_fp64)
          auto c = FIRFilterGPU<double>{};
          testDirac2(i, c);
        }
        {
          auto c = FIRFilterGPUAsync<float>{};
          testDirac2(i, c);
        }
        {
          // TODO test this on a machine that has support for doubles on the gpu (cl_khr_fp64)
          auto c = FIRFilterGPUAsync<double>{};
          testDirac2(i, c);
        }
        {
          auto c = FIRFilterGPUAsyncN<float>{};
          applySetup(c,{10});
          testDirac2(i, c);
        }
        {
          // TODO test this on a machine that has support for doubles on the gpu (cl_khr_fp64)
          auto c = FIRFilterGPUAsyncN<double>{};
          applySetup(c,{10});
          testDirac2(i, c);
        }
         */
        {
          auto c = PartitionnedFIRFilterGPUAsyncN<float>{};
          applySetup(c,{10});
          testDirac2(i, c);
        }
        {
          // TODO test this on a machine that has support for doubles on the gpu (cl_khr_fp64)
          auto c = PartitionnedFIRFilterGPUAsyncN<double>{};
          applySetup(c,{10});
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
          auto c = SubSampled<LatencySemantic::FirstNonZero, FFTConvolution<double, Tag>>{};
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
        testDiracPartitionned<float, Tag>(i);
        testDiracPartitionned<double, Tag>(i);
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
  
  auto c = mkRealTimeConvolution<double>();
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


