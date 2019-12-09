
namespace imajuscule
{

  template <typename Algo>
  struct Delayed {
    using T = typename Algo::FPT;
    using FPT = T;
    static constexpr int nCoefficientsFadeIn = Algo::nCoefficientsFadeIn;
    
    struct SetupParam : public Cost {
        using InnerParams = typename Algo::SetupParam;
        SetupParam(int d, InnerParams i)
        : delay(d)
        , innerParams(i)
        {}
        
        int delay;
        InnerParams innerParams;
        void logSubReport(std::ostream & os) const override {
            os << "Delayed by " << delay << std::endl;
            {
                IndentingOStreambuf i(os);
                innerParams.logSubReport(os);
            }
        }
    };
      void logComputeState(std::ostream & os) const {
          os << "Delayed by " << ring.size() << std::endl;
          IndentingOStreambuf i(os);
          algo.logComputeState(os);
      }

    void setup(SetupParam const & p) {
      ring.reset();
      ring.resize(p.delay);
        
      algo.setup(p.innerParams);
    }
    
    auto getEpsilon() const {
      return algo.getEpsilon();
    }
    
    auto isValid() const {
      return !ring.empty() && algo.isValid();
    }

    auto isZero() const {
      return algo.isZero();
    }

    void reset() {
      algo.reset();
      ring.resize(0);
    }
    void flushToSilence() {
      algo.flushToSilence();
      ring.reset();
    }
    
    auto getLatency() const {
      return ring.size() + algo.getLatency();
    }
    
    void setCoefficients(a64::vector<T> coeffs) {
      ring.reset();
      algo.setCoefficients(std::move(coeffs));
    }
    
    FPT step(FPT input) {
      auto oldest = *ring.cycleEnd();
      ring.feed(input);
      return algo.step(oldest);
    }
    template<typename FPT2>
    void stepAddVectorized(FPT2 const * const input_buffer,
                           FPT2 * output_buffer,
                           int nSamples)
    {
        for(int i=0; i<nSamples; ++i) {
            output_buffer[i] += step(input_buffer[i]);
        }
    }

    auto & getInner() {
      return algo;
    }
    
  private:
    cyclic<FPT> ring;
    Algo algo;
  };
}
