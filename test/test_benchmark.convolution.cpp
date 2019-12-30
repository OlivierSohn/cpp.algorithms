
template<typename Inputs, typename Conv>
auto test1( Inputs const & is, Conv & conv ) {
  typename Conv::FPT s = {};
  for(auto i : is) {
    s += conv.step(i);
  }
  return s;
}

template<typename Inputs, typename Conv>
auto test2( Inputs const & is, Conv & conv ) {
  std::vector<typename Conv::FPT> v;
  v.reserve(is.size());
  for(auto i : is) {
    v.push_back(conv.step(i));
  }
  return std::move(v);
}

namespace imajuscule {
template<typename T, typename Tag = fft::Fastest>
auto mkBruteThenScale(int firstSz, int nCoeffs) {
    using Scale = CustomScaleConvolution<FFTConvolutionIntermediate < PartitionnedFFTConvolutionCRTP<T, Tag> >>;
    using C = SplitConvolution<FIRFilter<T>,Scale>;
    using ScalingParam = typename Scale::SetupParam::ScalingParam;

    std::map<std::vector<Scaling>, double> results;
    
    {
        int count = 0;
        ScalingsIterator{
            firstSz,
            nCoeffs
        }.forEachScaling([nCoeffs, &count, &results](auto const & v){
            results.emplace(v, virtualCostPerSample(mkSimulation<Scale>(v, nCoeffs)));
            ++count;
        });
        if(results.size() != count) {
            throw std::logic_error("duplicate results");
        }
    }
    auto byCost = orderByValue(results);

    // write in the console the description of the latehandler
    {
        int const nDisplay = 1;
        analyzeScalings(firstSz, nCoeffs, byCost, nDisplay);
    }

    if(byCost.empty()) {
        throw std::logic_error("no optimal scaling");
    }
    C c;
    std::vector<Scaling> optimalScalings = byCost.front()->first;
    auto params = scalingsToParams<ScalingParam>(optimalScalings);
    c.setup(typename C::SetupParam{{},{params}});
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
  
  std::vector<std::vector<FPT>> dtDropped;
  static constexpr auto nDroppedMax = 12;
  dtDropped.resize(nDroppedMax);
  
  std::vector<int> coeffsCount;
  
  for(int powCoeffs=nDroppedMax; powCoeffs<nDroppedMax+4; ++powCoeffs) {
    coeffsCount.push_back(pow2(powCoeffs+1)-1);
  }
  
  FPT res{};

  for(auto nCoeffs : coeffsCount) {
    LG(INFO,"%d coefficients", nCoeffs);
    
    a64::vector<FPT> coeffs;
    coeffs.reserve(nCoeffs);
    for(int i=0; i<nCoeffs; ++i) {
      coeffs.push_back(cos(static_cast<FPT>(i)/100.f));
    }

    for(int firstSz = 1, nDropped = 0;
        nDropped < nDroppedMax;
        firstSz *= 2, ++nDropped)
    {
      auto c = mkBruteThenScale<FPT>(firstSz, coeffs.size());
      c.setCoefficients(coeffs);
      dtDropped[nDropped].push_back(measure_thread_cpu_one([&inputs, &c, &res](){
        res += test1(inputs,c);
      }).count());
    }

  }
  for(int i = 0; i<nDroppedMax; ++i) {
    LG(INFO,"%d dropped", i);
    auto plot = StringPlot(20,dtDropped[i].size());
    plot.drawLog(dtDropped[i], '+');
    plot.log();
  }
  
  for(int i = 0; i < coeffsCount.size(); ++i) {
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
    LG(INFO, "count coeffs: %d : min is dropped %d (%f)", coeffsCount[i], mi, m);
  }
  
  std::vector<std::vector<FPT>> dtDroppedTransposed;
  dtDroppedTransposed.resize(dtDropped[0].size()); // assuming they all have the same size
  for(int j=0; j<dtDropped[0].size(); ++j) {
    dtDroppedTransposed[j].resize(dtDropped.size());
    for(int i=0; i<dtDropped.size(); ++i) {
      dtDroppedTransposed[j][i] = dtDropped[i][j];
    }
  }

  for(int i = 0; i<coeffsCount.size(); ++i) {
    LG(INFO,"count coeffs: %d", coeffsCount[i]);
    auto plot = StringPlot(20,dtDroppedTransposed[i].size());
    plot.drawLog(dtDroppedTransposed[i], '+');
    plot.log();
  }

    std::cout << "dummy " << res << std::endl;
}

/*
 - Tuning of the number of dropped convolutions in
     SplitConvolution<
       FIRFilter<T>,
       CustomScaleConvolution<FFTConvolutionIntermediate < PartitionnedFFTConvolutionCRTP<T, Tag> >>
     >;
 
   The results show that on my platform, we should drop 5 scaled convolutions (sizes 1,2,4,8,16)
   and replace them by a single brute force convolution (of size 31).
 */
TEST(BenchmarkConvolutions, scaled_vs_brute) {
  test<double>();
}
