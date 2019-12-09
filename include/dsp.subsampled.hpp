
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
    static constexpr int nComputePhaseable = 0; // SubSampled phases will be set manually
    static constexpr int nCoefficientsFadeIn = scaleFadeSz::inSmallerUnits;
    static constexpr bool has_subsampling = true;

    using SetupParam = typename Algo::SetupParam;
 void logComputeState(std::ostream & os) const {
     os << "2x subsampled" << std::endl;
     IndentingOStreambuf i(os);
     algo.logComputeState(os);
 }

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
      
    std::array<int, 0> getComputePeriodicities() const {
        return {};
    }
    // in [0, getComputePeriodicity())
    std::array<int, 0> getComputeProgresses() const {
        return {};
    }
      void setComputeProgresses(std::array<int, 0> const & progresses) {
      }
      
    SubSampled() {
      resetStates();
    }
      
      void setup(SetupParam const & p) {
          algo.setup(p);
      }

    void setCoefficients(a64::vector<T> const & coeffs) {
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
        a64::vector<T> resampled;
        using audio::resampleSincBuffer;
        resampleSincBuffer(coeffs, 1, 2.0, resampled);
        assert(coeffs.size() == resampled.size() * 2);
      algo.setCoefficients(resampled);
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
    template<typename FPT2>
    void stepAddVectorized(FPT2 const * const input_buffer,
                           FPT2 * output_buffer,
                           int nSamples)
    {
        if(unlikely(isZero())) {
          return;
        }
        for(int i=0; i<nSamples; ++i) {
            output_buffer[i] += step(input_buffer[i]);
        }
    }
    template<typename FPT2>
    void stepAddInputZeroVectorized(FPT2 * output_buffer,
                                    int nSamples)
    {
        if(unlikely(isZero())) {
          return;
        }
        for(int i=0; i<nSamples; ++i) {
            output_buffer[i] += step({});
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
