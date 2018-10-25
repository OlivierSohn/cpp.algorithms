
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
      static void adaptCoefficients(a64::vector<F> & v) {
      }
    };
    
    template<LatencySemantic L, typename T>
    struct ConvolutionTraits<SubSampled<L, T>> {
      static constexpr bool supportsOddCountOfCoefficients = false;
      static constexpr bool evenIndexesAreApproximated = true;
      
      template <typename F>
      static void adaptCoefficients(a64::vector<F> & v) {
        averageNeighbours(v.begin(),v.end());
      }
    };
    
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
    void test(Convolution & conv, Coeffs const & coefficients, std::vector<bool> const & coeffMask) {
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
      std::vector<T> expectZero;
      std::vector<T> results;
      
      using Tr = ConvolutionTraits<Convolution>;
      
      if(conv.getLatency()) {
        expectZero.push_back(conv.step(1));
      }
      else {
        results.push_back(conv.step(1));
      }
      
      for(int i=1; i<conv.getLatency(); ++i) {
        expectZero.push_back(conv.step(0));
      }
      
      auto eps = conv.getEpsilon();
      for(auto j=0; j<expectZero.size(); ++j) {
        ASSERT_NEAR(0, expectZero[j], eps);
      }
      
      while(results.size() != coefficients.size()) {
        results.push_back(conv.step(0));
      }
      
      assert(results.size() == coefficients.size());

      for(auto j=0; j<results.size(); ++j) {
        if(!coeffMask[j]) {
          continue;
        }
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
    
    template<typename T, typename Tag>
    void testDiracFinegrainedPartitionned(int coeffs_index) {
      
      auto f = [](int part_size, auto const & coefficients)
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
          std::vector<bool> mask;
          mask.resize(coefficients.size(), true);
          test(conv, coefficients, mask);
        }
      };
      
      testPartitionned<T>(coeffs_index, f);
    }
    
    template<typename T, typename Tag>
    void testDiracPartitionned(int coeffs_index) {
      
      auto f = [](int part_size, auto const & coefficients){
        PartitionnedFFTConvolution<T,Tag> conv;
        
        applySetup(conv, {1.f,{part_size}});
        if(!conv.isValid()) {
          return;
        }
        conv.setCoefficients(coefficients);

        std::vector<bool> mask;
        mask.resize(coefficients.size(), true);
        test(conv, coefficients, mask);
      };
      
      testPartitionned<T>(coeffs_index, f);
    }
    
    
    template<typename Convolution>
    void testDirac2(int coeffs_index, Convolution & conv) {
      
      if(!conv.isValid()) {
        return;
      }
      
      using T = typename Convolution::FPT;
      
      auto coefficients = makeCoefficients<T>(coeffs_index);
      
      using Tr = ConvolutionTraits<Convolution>;
      if(!Tr::supportsOddCountOfCoefficients
         && 1 == coefficients.size() % 2) {
        return;
      }
      conv.setCoefficients(coefficients);
      
      Tr::adaptCoefficients(coefficients);

      std::vector<bool> mask;
      mask.reserve(coefficients.size());
      for(int i=0; i<coefficients.size(); ++i) {
        bool activate = true;
        if(Tr::evenIndexesAreApproximated && (0 == i%2)) {
          activate = false;
        }
        mask.push_back(activate);
      }

      test(conv, coefficients, mask);
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

      // crawl (A and B have non-overlapping output for a dirac)
      {
        using C = SplitConvolution<
        FIRFilter<double>,                      // A
        SubSampled<LatencySemantic::FirstNonZero, FFTConvolution<double, Tag>> // B
        >;
        
        auto c = C{};
        
        auto coefficients = makeCoefficients<double>(11);
        assert(coefficients.size() == 4 + 2*3);
        
        c.setSplit(4);
        
        using Tr = ConvolutionTraits<C>;
        c.setCoefficients(coefficients);
        
        auto it = coefficients.begin()+4;
        auto end = coefficients.end();
        averageNeighbours(it,end);
        
        std::vector<bool> mask;
        mask.reserve(coefficients.size());
        for(int j=0; j<coefficients.size(); ++j) {
          bool activate = true;
          if(j>=4 && (0 == (j-4)%2)) {
            activate = false;
          }
          mask.push_back(activate);
        }
        
        test(c, coefficients, mask);
      }
      
      // walk (A and B have 1 overlapping sample for a dirac input)
      {
        using C = SplitConvolution<
        FinegrainedPartitionnedFFTConvolution<double, Tag>,   // A
        SubSampled<
          LatencySemantic::DiracPeak,
          FinegrainedPartitionnedFFTConvolution<double, Tag>> // B
        >;
        
        // A has a latency of 7  (peak latency of 7, too)
        // B has a latency of 14 (peak latency of 15)
        // Hence the split should be at 15 - 7, and there is one overlapping sample.
        
        auto c = C{};
        applySetup(c,typename C::SetupParam{FinegrainedSetupParam{4,1,0},FinegrainedSetupParam{4,1,0}});

        auto coefficients = makeCoefficients<double>(11);
        assert(coefficients.size() == 4 + 2*3);
        
        c.setSplit(8); // B peak latency - A peak latency
        
        // so if we want to use a larger split, for example to have a better resolution for a longer time,
        // we should manually delay the inputs of B.
        
        using Tr = ConvolutionTraits<C>;
        c.setCoefficients(coefficients);
        
        auto it = coefficients.begin()+8;
        auto end = coefficients.end();
        averageNeighbours(it,end);
        
        std::vector<bool> mask;
        mask.reserve(coefficients.size());
        for(int j=0; j<coefficients.size(); ++j) {
          bool activate = true;
          if(j>=7 && (0 == (j-7)%2)) {
            activate = false;
          }
          mask.push_back(activate);
        }
        
        test(c, coefficients, mask);
      }
      
      // run (with delay)
      {
        using C = SplitConvolution<
        FinegrainedPartitionnedFFTConvolution<double, Tag>,            // A
        SubSampled<
          LatencySemantic::DiracPeak,
          Delayed<FinegrainedPartitionnedFFTConvolution<double, Tag>>> // B
        >;
        
        // A has a latency of 7  (peak latency of 7, too)
        // B has a latency of 2*delay + 14 (peak latency of 2*delay + 15)
        // Hence the split should be at 2*delay + 15 - 7, and there is one overlapping sample.
        // 2*delay = split + 7 - 15
        
        auto c = C{};
        // 4 is the min partition size to be valid
        applySetup(c,typename C::SetupParam
        {
          FinegrainedSetupParam{4,4,0},
          {
            1, // delay in downsampled units
            FinegrainedSetupParam{4,4,0}
          }
        });
        
        auto coefficients = makeCoefficients<double>(13);
        assert(coefficients.size() == 12);
        
        auto constexpr split = 10;
        
        c.setSplit(split);
        
        // so if we want to use a larger split, for example to have a better resolution for a longer time,
        // we should manually delay the inputs of B.
        
        using Tr = ConvolutionTraits<C>;
        c.setCoefficients(coefficients);
        
        auto it = coefficients.begin()+split;
        auto end = coefficients.end();
        averageNeighbours(it,end);
        
        std::vector<bool> mask;
        mask.reserve(coefficients.size());
        for(int j=0; j<coefficients.size(); ++j) {
          bool activate = true;
          if(j>=(split-1) && (0 == (j-(split-1))%2)) {
            activate = false;
          }
          mask.push_back(activate);
        }
        
        test(c, coefficients, mask);
      }
      
      
      
      for(int i=0; i<end_index; ++i) {
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
