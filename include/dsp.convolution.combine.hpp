
namespace imajuscule
{
  /*
   * Creates a convolution scheme by combining 2 convolution schemes,
   * where A handles ** early ** coefficients, and B handles ** late ** coefficients.
   */
  template<typename A, typename B>
  struct CombineConvolutions {
    using FPT = typename A::FPT;
    static_assert(std::is_same_v<FPT,typename B::FPT>);
    
    bool empty() const {
      return a.empty() && b.empty();
    }
    void clear() {
      a.clear();
      b.clear();
    }
    
    using SetupParamB = typename B::SetupParam;
    struct SetupParam {
      int split;
      SetupParamB setupB;
    };
    
    void applySetup(SetupParam const & p) {
      split = p.split;
      b.applySetup(p.setupB);
    }

    void setCoefficients(a64::vector<FPT> coeffs_) {
      auto [rangeA,rangeB] = splitAt(split, coeffs_);
      a.setCoefficients(rangeA.materialize());
      b.setCoefficients(rangeB.materialize());
      Assert(b.getLatency() == split + getLatency());
    }
    
    bool isValid() const {
      return a.isValid() && b.isValid();
    }

    void step(FPT val) {
      a.step(val);
      b.step(val);
    }
    
    auto get() const { return a.get() + b.get(); }

    auto getEpsilon() const {
      return a.getEpsilon() + b.getEpsilon();
    }

    auto getLatency() const {
      return a.getLatency();
    }
    
    auto & editA() {
      return a;
    }
    auto & editB() {
      return b;
    }
  private:
    int split; // we have 'split' ** early ** coefficients, the rest are ** late ** coefficients.
    A a;
    B b;
  };

  
  template<typename T>
  using NaiveConvolution = FIRFilter<T>;
  
  template<typename T, typename FFTTag = fft::Fastest>
  using RealTimeConvolution = CombineConvolutions<NaiveConvolution<T>,FinegrainedPartitionnedFFTConvolution<T, FFTTag>>;

  template<typename T, typename FFTTag = fft::Fastest>
  auto mkRealTimeConvolution(int split) {
    auto c = RealTimeConvolution<T, FFTTag>{};
    auto & b = c.editB();
    b.set_partition_size(split/2);
    c.applySetup({split,FinegrainedSetupParam{1,0}});
    return c;
  }

}
