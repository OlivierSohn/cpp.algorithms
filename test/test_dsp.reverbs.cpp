

static inline std::vector<double> mkDirac(int sz, double amplitude = 1.) {
  std::vector<double> res;
  res.resize(sz);
  if(res.empty()) {
    throw std::logic_error("0-size");
  }
  res[0] = amplitude;
  return res;
}

static inline std::vector<double> mkCoefficientsRamp(int sz) {
  std::vector<double> res;
  res.reserve(sz);
  for(int i=0; i<sz; ++i) {
    res.push_back(1.-(i/(double)sz));
  }
  return res;
}

// the coefficients should tend to 0, to mimic real responses
// (else, with scaling, the end of the dirac's response will be wrong)
static inline std::vector<double> mkCoefficientsTriangle(int sz) {
  std::vector<double> res;
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

static inline void scaleVec(double coeff, std::vector<double> & v) {
  std::for_each(v.begin(), v.end(), [coeff](auto & val) { val *= coeff; });
}

TEST(Reverbs, dirac) {
  using namespace imajuscule;

  Reverbs<1> rs;
  

  std::vector<int> audio_cb_sizes{1,2,3,4,7,8,9,63,64,65,256,1024,4096};
  
  for(auto cb_sz:audio_cb_sizes)
  {
    LG(INFO, "cb_sz = %d", cb_sz);
    // by default, reverb is inactive.
    {
      std::vector<double> const input = mkDirac(4);
      auto inputCopy = input;
      rs.apply(inputCopy.data(), inputCopy.size());
      ASSERT_EQ(inputCopy, input);
    }
    
    // for 0-length responses, reverb is inactive.
    {
      rs.setConvolutionReverbIR({}, 1, cb_sz, 44100.);
      std::vector<double> const input = mkDirac(4);
      auto inputCopy = input;
      rs.apply(inputCopy.data(), inputCopy.size());
      ASSERT_EQ(inputCopy, input);
    }
  
// TODO test larger sizes, up-to when all scales are used
    for(int sz = 1; sz < 1000; ++sz) {
      std::vector<double> const input = mkDirac(sz);
      auto inputCopy = input;
      
      auto coeffs = mkCoefficientsTriangle(sz);
      rs.setConvolutionReverbIR(coeffs, 1, cb_sz, 44100.);
      
      std::vector<double> diff1;
      diff1.reserve(inputCopy.size());

      rs.apply(inputCopy.data(), 1);
      auto scale = coeffs[0]/inputCopy[0];
      inputCopy[0] *= scale;
      diff1.push_back(inputCopy[0] - coeffs[0]);

      for(int i=1; i<inputCopy.size(); i++) {
        rs.apply(&inputCopy[i], 1);
        inputCopy[i] *= scale;
        
        diff1.push_back(std::abs(inputCopy[i] - coeffs[i]));
      }
      
      std::vector<std::pair<double, int>> diff;
      diff.reserve(diff1.size());
      for(int i=0; i<diff1.size(); ++i) {
        diff.emplace_back(diff1[i], i);
      }
      std::sort(diff.begin(), diff.end(), std::greater<>());
      int const n_scales = rs.countScales();
      
      ASSERT_LT(0, n_scales);
      
      //int const totalFades = (pow2(n_scales-1)-1) * scaleFadeSz::inSmallerUnits;
      
      switch(n_scales) {
        case 1:
          ASSERT_GT(0.00001, diff[0].first);
          break;
        case 2:
          ASSERT_GT(0.05, diff[0].first);
          break;
        case 3:
          ASSERT_GT(0.08, diff[0].first);
          break;
        case 4:
          ASSERT_GT(0.1, diff[0].first);
          break;
        default:
          ASSERT_TRUE(false);
          break;
      }
    }
    
    // reverb is inactive when disabled
    {
      rs.disable();
      std::vector<double> const input = mkDirac(4);
      auto inputCopy = input;
      rs.apply(inputCopy.data(), inputCopy.size());
      ASSERT_EQ(inputCopy, input);
    }
  }
}

