
namespace imajuscule
{
  /*
   * Creates a convolution scheme by combining 2 convolution schemes,
   * where A handles ** early ** coefficients, and B handles ** late ** coefficients.
   */
  template<typename A, typename B>
  struct CombineConvolutions {
    using FPT = typename A::FPT;
    using EarlyHandler = A;
    using LateHandler = B;

    static_assert(std::is_same_v<FPT,typename B::FPT>);

    bool empty() const {
      return a.empty() && b.empty();
    }
    void clear() {
      a.clear();
      b.clear();
    }

    using SetupParam = typename B::SetupParam;

    void set_partition_size(int sz) {
      b.set_partition_size(sz);
      split = b.getLatency() - a.getLatency();
    }

    void setCoefficients(a64::vector<FPT> coeffs_) {
      auto [rangeA,rangeB] = splitAt(split, coeffs_);
      a.setCoefficients(rangeA.materialize());
      b.setCoefficients(rangeB.materialize());
    }

    void applySetup(SetupParam const & p) {
      b.applySetup(p);
    }

    auto getGranularMinPeriod() const {
      return b.getGranularMinPeriod();
    }

    bool isValid() const {
      return a.isValid() && b.isValid();
    }

    void step(FPT val) {
      a.step(val);
      b.step(val);
    }

    auto get() const {
      return a.get() + b.get();
    }

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


  template<typename T>
  struct PartitionAlgo< RealTimeConvolution<T> > {
    using NonAtomicConvolution = RealTimeConvolution<T>;
    using Delegate = typename NonAtomicConvolution::LateHandler;

    using SetupParam = typename Delegate::SetupParam;
    using PartitionningSpec = PartitionningSpec<SetupParam>;
    using PartitionningSpecs = PartitionningSpecs<SetupParam>;

    static PartitionningSpecs run(int n_channels, int n_audio_frames_per_cb, int size_impulse_response) {
      return PartitionAlgo<Delegate>::run(n_channels, n_audio_frames_per_cb, size_impulse_response);
    }
  };
}
