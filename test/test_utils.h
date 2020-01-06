namespace imajuscule {

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

template<typename T>
std::unique_ptr<T> extractFromEnd(std::vector<std::unique_ptr<T>> &c, T*o) {
  for(auto i = c.rbegin(), e = c.rend(); i!=e; ++i) {
    auto & p = *i;
    if(p.get() != o) {
      continue;
    }
    std::unique_ptr<T> res;
    res.swap(p);
    c.erase(std::next(i).base());
    return res;
  }
  return {};
}

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

// the coefficients should tend to 0, to mimic real responses
// (else, with scaling, the end of the dirac's response will be wrong)
static inline a64::vector<double> mkCoefficientsTriangle(int sz) {
    a64::vector<double> res;
    res.reserve(sz);
    constexpr int period = 100;
    for(int i=0; i<sz; ++i) {
        auto j = i%(2*period);
        if(j < period) {
            res.push_back(1.-(j/(double)period));
        }
        else {
            res.push_back((j-period)/(double)period);
        }
    }
    
    constexpr int fadeout = 40;
    if(sz > fadeout) {
        for(int i=0; i<fadeout; ++i) {
            res[res.size()-1-i] *= (i+1) / static_cast<double>(fadeout);
        }
    }
    return res;
}

std::vector<Scaling> mkNaiveScaling(int firstSz, int const countCoeffs);
std::vector<Scaling> mkBetterScaling(int firstSz, int const countCoeffs);

template<typename T, typename FFTTag>
auto mkRealTimeConvolutionSubsampled(std::vector<Scaling> const & v, int partitionSize) {
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


template<typename T, typename FFTTag>
auto mkRealTimeConvolution(std::vector<Scaling> const & v, int partitionSize) {
  using C = ZeroLatencyScaledFineGrainedPartitionnedConvolution<T, FFTTag>;
  using ScalingParam = typename C::SetupParam::AParam::BParam::ScalingParam;
  
  auto scalingParams = scalingsToParams<ScalingParam>(v);
  auto c = C{};
  c.setup(typename C::SetupParam
  {
    {
        {},
        {scalingParams}
    },
      FinegrainedSetupParam{
          partitionSize, // partition size
          partitionSize*1000, // multiplication group size
          0 // phase
      }
    
  }
             );
  return c;
}


template<typename T, typename Tag>
auto mkRealTimeConvolution2(std::vector<Scaling> const & v,
                            int partitionSize,
                            int partitionCount) {
  using C = Convolution<
    AlgoSplitConvolution<
      AlgoSplitConvolution <
          AlgoFIRFilter<T, Tag>,
          AlgoCustomScaleConvolution<AlgoFFTConvolutionIntermediate < AlgoPartitionnedFFTConvolutionCRTP<T, Tag> >>
        >,
      AlgoFinegrainedPartitionnedFFTConvolution<T, Tag>
  >>;
  using ScalingParam = typename C::SetupParam::AParam::BParam::ScalingParam;
  
  auto scalingParams = scalingsToParams2<ScalingParam>(v);
  auto c = C{};
  c.setup(typename C::SetupParam
  {
    {
        {},
        {scalingParams}
    },
      FinegrainedSetupParam2{
          partitionSize, // partition size
          partitionCount,
          partitionSize*1000, // multiplication group size
          0 // phase
      }
    
  }
             );
  return c;
}
}
