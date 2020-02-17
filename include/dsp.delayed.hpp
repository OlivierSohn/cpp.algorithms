
namespace imajuscule
{

  template <typename Algo>
  struct Delayed {
    using T = typename Algo::FPT;
    using FPT = T;

    static constexpr int nCoefficientsFadeIn = Algo::nCoefficientsFadeIn;
    static constexpr bool step_can_error = Algo::step_can_error;

    struct SetupParam : public Cost {
        static constexpr int nCoefficientsFadeIn = Algo::nCoefficientsFadeIn;

        using InnerParams = typename Algo::SetupParam;
        SetupParam(int d, InnerParams i)
        : delay(d)
        , delayed(i)
        {}
        
        int delay;
        InnerParams delayed;
        
        void logSubReport(std::ostream & os) const override {
            os << "Delayed by " << delay << std::endl;
            {
                IndentingOStreambuf i(os);
                delayed.logSubReport(os);
            }
        }
        
        bool handlesCoefficients() const {
            return delayed.handlesCoefficients();
        }
        
        Latency getImpliedLatency() const {
          Assert(handlesCoefficients());
          return Latency(delay) + delayed.getImpliedLatency();
        }

        void adjustWork(int targetNCoeffs) {
            if(delayed.handlesCoefficients()) {
                delayed.adjustWork(targetNCoeffs);
            }
            if(!handlesCoefficients()) {
                setCost(0.);
            }
        }
    };
      void logComputeState(std::ostream & os) const {
          os << "Delayed by " << ring.size() << std::endl;
          {
              IndentingOStreambuf i(os);
              algo.logComputeState(os);
          }
      }

    void setup(SetupParam const & p) {
      ring.reset();
      ring.resize(p.delay);
        
      algo.setup(p.delayed);
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
    
      bool handlesCoefficients() const {
          return algo.handlesCoefficients();
      }
    Latency getLatency() const {
      Assert(handlesCoefficients());
      return Latency(ring.size()) + algo.getLatency();
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
    
      int const getBiggestScale() const {
          return 2*algo.getBiggestScale();
      }
      
  private:
    cyclic<FPT> ring;
    Algo algo;
  };
}
