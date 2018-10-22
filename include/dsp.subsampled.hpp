
namespace imajuscule
{
  enum class LatencySemantic {
    // the delay between an input dirac and the first non-zero output:
    FirstNonZero,
    // the delay between an input dirac and the peak output:
    DiracPeak
  };

  /* The impulse response is downsampled by a factor of 2. */
  template <LatencySemantic Lat, typename Algo>
  struct SubSampled {
    using T = typename Algo::FPT;
    using FPT = T;
    
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
      if constexpr (Lat == LatencySemantic::DiracPeak) {
        return 1 + 2 * algo.getLatency();
      }
      return 2 * algo.getLatency();
    }
    
    SubSampled() {
      reset();
    }
    
    void setCoefficients(a64::vector<T> coeffs) {
      reset();
      
      assert(coeffs.size() % 2 == 0);

      auto writeIt = coeffs.begin();
      {
        auto readIt = coeffs.begin();
        
        auto const end = coeffs.end();
        while(readIt+1 < end)
        {
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
      clock = true;
      prevInput = 0;
      prevOutput = 0;
      algo.reset();
    }

  private:
    Algo algo;
    bool clock; // not needed in branchless mode
    FPT prevOutput;
    FPT prevInput; // not needed in branchless mode

  };
  
}
