
namespace imajuscule
{

  template <typename Algo>
  struct Delayed {
    using T = typename Algo::FPT;
    using FPT = T;
    static constexpr int nCoefficientsFadeIn = Algo::nCoefficientsFadeIn;
    
    struct SetupParam {
      int delay;
      typename Algo::SetupParam innerParams;
    };
    
    // instead of calling that directly, use applySetup
    void setTheDelay(int d) {
      ring.reset();
      ring.resize(d);
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
    
    auto & getInner() {
      return algo;
    }
    
  private:
    cyclic<FPT> ring;
    Algo algo;
  };
}
