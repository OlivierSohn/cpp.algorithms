
namespace imajuscule
{
  enum class LatencySemantic {
    // the delay between an input dirac and the first non-zero output:
    FirstNonZero,
    // the delay between an input dirac and the peak output:
    DiracPeak
  };
  
  /*
   * To have a smooth transition between different sampling frequencies,
   * we cross-fade the impulse response coefficients.
   */
  struct scaleFadeSz {
    // expressed in number of periods at the higher frequency
    static int constexpr inSmallerUnits = 500;
    static_assert(inSmallerUnits % 2 == 0);
    // expressed in number of periods at the lower frequency
    static int constexpr inBiggerUnits = inSmallerUnits/2;
  };
  constexpr int nMaxScales = 4;

  // we want to avoid approximations:
  constexpr bool subSamplingAllowsEvenNumberOfCoefficients = false;

  /* The impulse response is downsampled by a factor of 2. */
  template <LatencySemantic Lat, typename Algo>
  struct SubSampled {
    using T = typename Algo::FPT;
    using FPT = T;
    static constexpr int nCoefficientsFadeIn = scaleFadeSz::inSmallerUnits;
    
    using SetupParam = typename Algo::SetupParam;
 
    auto getEpsilon() const {
      return algo.getEpsilon();
    }
    
    auto isValid() const {
      return algo.isValid();
    }
    
    auto isZero() const {
      return algo.isZero();
    }
    
    auto getLatency() const {
      int const res = 2 * algo.getLatency();
      if constexpr (Lat == LatencySemantic::DiracPeak) {
        return 1 + res;
      }
      return res;
    }
    
    SubSampled() {
      resetStates();
    }
    
    void setCoefficients(a64::vector<T> coeffs) {
      resetStates();

      if constexpr (subSamplingAllowsEvenNumberOfCoefficients) {
        if(coeffs.size()%2) {
          // duplicate last coefficient.
          coeffs.push_back(coeffs.back());
        }
      }
      else {
        assert(coeffs.size() % 2 == 0);
      }

      auto writeIt = coeffs.begin();
      {
        auto readIt = coeffs.begin();
        
        auto const end = coeffs.end();
        while(readIt+1 < end)
        {
            // UGLY!!! TODO use sinc interpolation
          *writeIt = 0.5 * (*readIt + *(readIt+1));
          
          ++writeIt;
          readIt += 2;
        }
        
        if(readIt != end) {
          throw std::logic_error("cannot subsample an odd count of coefficients");
        }
      }

      coeffs.resize(std::distance(coeffs.begin(), writeIt));

      algo.setCoefficients(std::move(coeffs));
    }
    
    FPT step(FPT input) {
      clock = !clock;
      
      // This branch will be badly predicted all the time
      // because it changes every time it is taken
      if(clock) {
        prevInput = input;
        return prevOutput;
      }
      else {
        // we don't half the input because the coefficients were halved in setCoefficients
        auto sumInputs = prevInput + input;
        auto output = algo.step(sumInputs);
        auto res = 0.5f * (output + prevOutput);
        prevOutput = output;
        return res;
      }
    }
    
    auto & getInner() {
      return algo;
    }

    void reset() {
      resetStates();
      algo.reset();
    }

    void flushToSilence() {
      resetStates();
      algo.flushToSilence();
    }

  private:
    Algo algo;
    bool clock;
    FPT prevOutput;
    FPT prevInput;
    
    void resetStates() {
      clock = true;
      prevInput = 0;
      prevOutput = 0;
    }
  };
  
}
