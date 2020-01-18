
namespace imajuscule {

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

template<LatencySemantic Lat, typename InnerParam>
struct SubSampledSetupParam : public Cost {
    static constexpr int nCoefficientsFadeIn = scaleFadeSz::inSmallerUnits;

    SubSampledSetupParam(InnerParam const & i)
    : subsampled(i)
    {}
    
    bool handlesCoefficients() const {
        return subsampled.handlesCoefficients();
    }
    Latency getImpliedLatency() const {
        Assert(handlesCoefficients());
        int const res = 2 * subsampled.getImpliedLatency().toInteger();
        if constexpr (Lat == LatencySemantic::DiracPeak) {
            return Latency(1 + res);
        }
        return Latency(res);
    }
    
    void adjustWork(int targetNCoeffs) {
        if(targetNCoeffs % 2) {
            ++targetNCoeffs;
        }
        
        subsampled.adjustWork(targetNCoeffs/2);
        if(!handlesCoefficients()) {
            setCost(0.);
        }
        else {
            setCost(subsampled.getCost()/2.);
        }
    }
    
    void logSubReport(std::ostream & os) const override
    {
        os << "2x subsampled" << std::endl;
        {
            IndentingOStreambuf i(os);
            subsampled.logSubReport(os);
        }
    }

    InnerParam subsampled;
};

  /* The impulse response is downsampled by a factor of 2. */
  template <LatencySemantic Lat, typename Algo>
  struct SubSampled {
    using T = typename Algo::FPT;
    using FPT = T;
    static constexpr int nComputePhaseable = 0; // SubSampled phases will be set manually
    static constexpr int nCoefficientsFadeIn = scaleFadeSz::inSmallerUnits;
    static constexpr bool has_subsampling = true;
    static constexpr bool step_can_error = Algo::step_can_error;

    using SetupParam = SubSampledSetupParam<Lat, typename Algo::SetupParam>;
 void logComputeState(std::ostream & os) const {
     os << "2x subsampled" << std::endl;
     {
         IndentingOStreambuf i(os);
         algo.logComputeState(os);
     }
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
    
      bool handlesCoefficients() const {
          return algo.handlesCoefficients();
      }
      
    Latency getLatency() const {
      Assert(handlesCoefficients());
      int const res = 2 * algo.getLatency().toInteger();
      if constexpr (Lat == LatencySemantic::DiracPeak) {
        return Latency(1 + res);
      }
      return Latency(res);
    }
      
    std::array<int, nComputePhaseable> getComputePeriodicities() const {
        return {};
    }
    // in [0, getComputePeriodicity())
    std::array<int, nComputePhaseable> getComputeProgresses() const {
        return {};
    }
      void setComputeProgresses(std::array<int, nComputePhaseable> const & progresses) {
      }
      
    SubSampled() {
      resetStates();
    }
      
      void setup(SetupParam const & p) {
          algo.setup(p.subsampled);
      }

    void setCoefficients(a64::vector<T> /* by copy */ coeffs) {
      resetStates();

        if(coeffs.size()%2) {
            coeffs.push_back(T{});
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
        auto res = 0.5f * (output + prevOutput); // TODO there is probably much better we can do here, like sinc interpolation, at the expense of bigger latency.
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
