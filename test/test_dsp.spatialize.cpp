
namespace imajuscule {
  namespace testspatialize {
    constexpr auto end_index = 4;
    
    template<typename T>
    a64::vector<T> makeCoefficients(int coeffs_index) {
      switch(coeffs_index) {
        case 0: return {{ .9,.8,.7,.6,.3,.2,.1,0. }};
        case 1: return {{ +1. }};
        case 2: return {{ -1. }};
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
      
      int lat = spatialize.getLatency();
      
      std::vector<std::array<T, nEars>> results;
      
      // feed a dirac
      results.emplace_back();
      std::fill(results.back().begin(), results.back().end(), 1);

      while(results.size() < lat + ear_signals[0].size()) {
        results.emplace_back();
        std::fill(results.back().begin(), results.back().end(), 0);
      }

      // process the dirac
      for(auto & r : results) {
        spatialize.step(r.data(), 0., 1.);
      }
      
      auto eps = spatialize.getEpsilon();
      for(int j=0; j<lat; ++j) {
        int i=0;
        for(auto r : results[j]) {
          ASSERT_NEAR(0, r, eps);
          ++i;
        }
      }
      for(int j=lat; j<results.size(); ++j) {
        int i=0;
        for(auto r : results[j]) {
          ASSERT_NEAR(ear_signals[i][j-lat], r, eps);
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
        
        std::pair<int,typename Convolution::SetupParam> setup {part_size,{}};
        
        {
          constexpr auto nOutMono = 1;
          audio::Spatializer<nOutMono, Convolution> spatialized;
          
          spatialized.addSourceLocation({{coefficients}}, setup);
          
          test(spatialized, {{coefficients}});
        }
        {
          constexpr auto nOutStereo = 2;
          audio::Spatializer<nOutStereo, Convolution> spatialized;
          
          // no cross-talk :
          // source1 has only left component
          // source2 has only right component
          spatialized.addSourceLocation({{coefficients, zero_coeffs}}, setup);
          spatialized.addSourceLocation({{zero_coeffs, coefficients}}, setup);
          
          test(spatialized, {{coefficients, coefficients}});
        }
        
        {
          constexpr auto nOutStereo = 2;
          audio::Spatializer<nOutStereo, Convolution> spatialized;
          
          opposite_coeffs = coefficients;
          for(auto & o : opposite_coeffs) {
            o = -o;
          }
          // extreme cross-talk :
          // source1 right component is the opposite of source 2
          // source2 right component is the opposite of source 1
          spatialized.addSourceLocation({{coefficients, opposite_coeffs}}, setup);
          spatialized.addSourceLocation({{opposite_coeffs, coefficients}}, setup);
          
          test(spatialized, {{zero_coeffs,zero_coeffs}});
        }
        {
          constexpr auto nOutStereo = 2;
          audio::Spatializer<nOutStereo, Convolution> spatialized;
          
          opposite_coeffs = coefficients;
          for(auto & o : opposite_coeffs) {
            o = -o;
          }
          
          spatialized.addSourceLocation({{coefficients, opposite_coeffs}}, setup);
          spatialized.addSourceLocation({{zero_coeffs, coefficients}}, setup);
          
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

