

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
    for(int sz = 1; sz < 10000; ++sz) {
      std::vector<double> const input = mkDirac(sz);
      auto inputCopy = input;
      
      auto coeffs = mkCoefficientsRamp(sz);
      std::reverse(coeffs.begin(), coeffs.end());
      for(auto & v:coeffs) {
        v+=1;
      }
      rs.setConvolutionReverbIR(coeffs, 1, cb_sz, 44100.);
      
      std::vector<double> diff;
      diff.reserve(inputCopy.size());

      rs.apply(inputCopy.data(), 1);
      auto scale = coeffs[0]/inputCopy[0];
      inputCopy[0] *= scale;
      diff.push_back(inputCopy[0] - coeffs[0]);

      for(int i=1; i<inputCopy.size(); i++) {
        rs.apply(&inputCopy[i], 1);
        inputCopy[i] *= scale;
        
        diff.push_back(inputCopy[i] - coeffs[i]);
        // uncomment to debug sample by sample:
        //*
      }
      for(int i=0; i<inputCopy.size(); i++)
      {
        //*/
        if(std::abs(inputCopy[i] - coeffs[i]) / (std::abs(inputCopy[i]) + std::abs(coeffs[i])) > 0.0001) {
          LG(ERR, "size %d at index %d: %f %f", sz, i, inputCopy[i]/scale, coeffs[i]/scale);
        }
        ASSERT_FLOAT_EQ(coeffs[i],inputCopy[i]);
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

