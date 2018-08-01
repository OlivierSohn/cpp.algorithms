
template<typename Inputs, typename Conv>
auto test1( Inputs const & is, Conv & conv ) {
  typename Conv::FPT s = {};
  for(auto i : is) {
    conv.step(i);
    s += conv.get();
  }
  return s;
}

template<typename Inputs, typename Conv>
auto test2( Inputs const & is, Conv & conv ) {
  std::vector<typename Conv::FPT> v;
  v.reserve(is.size());
  for(auto i : is) {
    conv.step(i);
    v.push_back(conv.get());
  }
  return std::move(v);
}

namespace imajuscule {
  template<typename FPT, typename FFTTag = fft::Fastest>
  auto mkBruteThenScale(int nDropped) {
    using C = SplitConvolution<FIRFilter<FPT>,ScaleConvolution<FFTConvolutionCore<FPT, FFTTag>>> ;
    C c;
    applySetup(c, typename C::SetupParam{{},{nDropped}});
    return c;
  }
}

template<typename FPT>
void test() {
  using namespace imajuscule;
  using namespace imajuscule::profiling;
  using namespace std::chrono;
  
  std::vector<FPT> inputs;
  constexpr auto inputSz = 1000000;
  inputs.reserve(inputSz);
  for(int i=0; i<inputSz; ++i) {
    inputs.push_back(sin(i));
  }
  
  std::vector<FPT> dtBrute, epsilonBrute;
  std::vector<FPT> dtScale, epsilonScale;
  std::vector<std::vector<FPT>> dtDropped;
  static constexpr auto nDroppedMax = 10;
  dtDropped.resize(nDroppedMax);
  
  std::vector<int> coeffs;
  
  for(int powCoeffs=0; powCoeffs<16; ++powCoeffs) {
    /*
    for(int i=0; i<2; ++i) {
      auto nCoeffs = pow2(powCoeffs);
      if(i) {
        nCoeffs *= 1.414;
      }
      coeffs.push_back(nCoeffs);
    }*/
    // count of cefficients used when doing 0-latency split convolutions:
    coeffs.push_back(pow2(powCoeffs+1)-1);
    //coeffs.push_back(pow2(powCoeffs));
  }
  /*
  coeffs.resize(511-255);
  std::iota(coeffs.begin(), coeffs.end(), 255);
  */
  /*
   coeffs.resize(128);
   std::iota(coeffs.begin(), coeffs.end(), 1);
*/
  
  for(auto nCoeffs : coeffs) {
    LG(INFO,"%d coefficients", nCoeffs);
    
    a64::vector<FPT> coeffs;
    for(int i=0; i<nCoeffs; ++i) {
      coeffs.push_back(cos(static_cast<FPT>(i)/100.f));
    }
    coeffs.reserve(nCoeffs);

    for(int nDropped = 0; nDropped < nDroppedMax; ++nDropped) {
      auto c = mkBruteThenScale<FPT>(nDropped);
      c.setCoefficients(coeffs);
      FPT res;
      dtDropped[nDropped].push_back(measure_one<high_resolution_clock>([&inputs, &c, &res](){
        res = test1(inputs,c);
      }));
      
      std::cout << res << std::endl;
    }

    /*
    {
      //ScaleConvolution<FFTConvolutionCore<FPT>> c(0);
      FFTConvolution<FPT> c;
      c.setCoefficients(coeffs);
      epsilonScale.push_back(c.getEpsilon());
      
      FPT res;
      dtScale.push_back(measure_one<high_resolution_clock>([&inputs, &c, &res](){
        res = test1(inputs,c);
      }));
      
      std::cout << res << std::endl;
    }
    {
      FIRFilter<FPT> c;
      c.setCoefficients(coeffs);
      epsilonBrute.push_back(c.getEpsilon());

      FPT res;
      dtBrute.push_back(measure_one<high_resolution_clock>([&inputs, &c, &res](){
        res = test1(inputs,c);
      }));
      
      std::cout << res << std::endl;
    }*/
  }
  
  {
    auto plot = StringPlot(20,dtBrute.size());
    plot.drawLog(dtBrute, '*');
    plot.drawLog(dtScale);
    plot.log();
  }
  {
    auto plot = StringPlot(20,dtBrute.size());
    plot.drawLog(dtScale);
    plot.log();
  }
  {
    auto plot = StringPlot(20,dtBrute.size());
    plot.draw(dtScale);
    plot.log();
  }

  for(int i=0; i<dtBrute.size(); ++i) {
    std::cout << coeffs[i] << " : " << dtBrute[i]/1e9 << " " << dtScale[i]/1e9 << std::endl;
  }
  {
    auto plot = StringPlot(20,dtBrute.size());
    plot.drawLog(epsilonBrute, '*');
    plot.drawLog(epsilonScale);
    plot.log();
  }
  for(int i=0; i<dtBrute.size(); ++i) {
    std::cout << coeffs[i] << " : " << epsilonBrute[i] << " " << epsilonScale[i] << std::endl;
  }
  
  for(int i = 0; i<nDroppedMax; ++i) {
    LG(INFO,"%d dropped", i);
    auto plot = StringPlot(20,dtDropped[i].size());
    plot.drawLog(dtDropped[i], '+');
    plot.log();
  }
  
  for(int i = 0; i < coeffs.size(); ++i) {
    auto m = std::numeric_limits<FPT>::max();
    auto mi = -1;
    auto j = -1;
    for(auto const & d : dtDropped) {
      ++j;
      if(m < d[i]) {
        continue;
      }
      m = d[i];
      mi = j;
    }
    LG(INFO, "coeff %d : min is dropped %d (%f)", coeffs[i], mi, m);
  }
  
  std::vector<std::vector<FPT>> dtDroppedTransposed;
  dtDroppedTransposed.resize(dtDropped[0].size()); // assuming they all have the same size
  for(int j=0; j<dtDropped[0].size(); ++j) {
    dtDroppedTransposed[j].resize(dtDropped.size());
    for(int i=0; i<dtDropped.size(); ++i) {
      dtDroppedTransposed[j][i] = dtDropped[i][j];
    }
  }

  for(int i = 0; i<coeffs.size(); ++i) {
    LG(INFO,"coeff %d", coeffs[i]);
    auto plot = StringPlot(20,dtDroppedTransposed[i].size());
    plot.drawLog(dtDroppedTransposed[i], '+');
    plot.log();
  }
}

/*
 - Tuning of the number of dropped convolutions in
     SplitConvolution<
       FIRFilter<FPT>,
       ScaleConvolution<FFTConvolutionCore<FPT, FFTTag>>
     >;
 
   The results show that on my platform, we should drop 5 scaled convolutions (sizes 1,2,4,8,16)
   and replace them by a single brute force convolution (of size 31).
 
 - other tests, commented out.
 */
TEST(BenchmarkConvolutions, scaled_vs_brute) {
  test<double>();
}
