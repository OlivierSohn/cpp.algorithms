
namespace imajuscule
{
  /*
   * Creates a convolution scheme by combining 2 convolution schemes,
   * where A handles ** early ** coefficients, and B handles ** late ** coefficients.
   */
  template<typename A, typename B>
  struct SplitConvolution {
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

  private:
    int split; // we have 'split' ** early ** coefficients, the rest are ** late ** coefficients.
    A a;
    B b;
  };


  /*
   * Creates a 0-latency convolution scheme by combining 'n' sub-schemes
   * of latencies 2^k-1, where k is in [0,n-1].
   */
  template<typename A>
  struct ScaleConvolution {
    using FPT = typename A::FPT;

    bool empty() const {
      return v.empty() || std::all_of(v.begin(), v.end(), [](auto & e) -> bool { return e.empty(); });
    }
    void clear() {
      v.clear();
      delayBuffer.resize(0);
    }

    /*
     *   Unless the count of coefficients is of the form
     *     2^n - 1, some padding occurs.
     */
    void setCoefficients(a64::vector<FPT> coeffs_) {
      clear();
      auto s = coeffs_.size();
      if( s == 0 ) {
        return;
      }
      auto n = power_of_two_exponent(ceil_power_of_two(s+1));
      v.resize(n);
      auto it = coeffs_.begin();
      auto end = coeffs_.end();
      for(int i=0; i<n; ++i) {
        Assert(it <= end);
        auto start = it;
        auto sizeBlock = pow2(i);
        it += sizeBlock;
        if(it > end) {
          Assert(i == n-1);
          auto withPadding = a64::vector<FPT>{start,end};
          // We could do naive padding, up-to sizeBlock, but to be
          // more efficient, memory- and time-wise, we pad
          // up to the power of 2 following std::distance(start,end),
          // and emulate the missing latency by buffering the values.
          auto nextPowOf2 = ceil_power_of_two(withPadding.size());
          withPadding.resize(nextPowOf2);
          Assert(sizeBlock >= nextPowOf2);
          delayBuffer.resize(sizeBlock-nextPowOf2);
          LG(INFO,"ScaleConvolution : delay buffer of period %d for the last convolution", delayBuffer.period());
          v[i].setCoefficients(withPadding);
          sizeBlock = withPadding.size();
        }
        else {
          v[i].setCoefficients({start,it});
        }
        if(v[i].getLatency() != sizeBlock - 1) {
          // This breaks the class logic, and would lead to wrong results.
          //   (for example, ScaleConvolution<FinegrainedPartitionnedFFTConvolution<float>> is not usable with this class.)
          LG(ERR,"ScaleConvolution is applied to a type that doesn't respect the latency constraints.");
          delayBuffer.resize(0);
          return;
        }
      }
    }

    bool isValid() const {
      return std::all_of(v.begin(), v.end(), [](auto & e) -> bool { return e.isValid(); });
    }

    void step(FPT val) {
      if(delayBuffer.zeroPeriod()) {
        // step all
        std::for_each(v.begin(), v.end(), [val](auto & e) { e.step(val); });
      }
      else if(!v.empty()) {
        // step all but the last one
        std::for_each(v.begin(), v.end()-1, [val](auto & e) { e.step(val); });
        // step the last one using the delayed value
        v.back().step(*delayBuffer.cycleEnd());
        delayBuffer.feed(val);
      }
    }

    auto get() const {
      return std::transform_reduce(v.begin(), v.end(), FPT{}, std::plus<>(), [](auto & e) { return e.get(); });
    }

    auto getEpsilon() const {
      return std::transform_reduce(v.begin(), v.end(), FPT{}, std::plus<>(), [](auto & e) { return e.getEpsilon(); });
    }

    auto getLatency() const { return 0; }

  private:
    std::vector<A> v;
    cyclic<FPT> delayBuffer; // used for the last convolution, to optimize padding.
                             // It emulates the effect of a bigger latency when using a bigger padding.
  };


  template<typename T>
  using ZeroLatencyBruteConvolution = FIRFilter<T>;

  template<typename T, typename FFTTag = fft::Fastest>
  using ZeroLatencyBruteFineGrainedPartitionnedConvolution = SplitConvolution<ZeroLatencyBruteConvolution<T>,FinegrainedPartitionnedFFTConvolution<T, FFTTag>>;

  template<typename T, typename FFTTag = fft::Fastest>
  using ZeroLatencyScaledFineGrainedPartitionnedConvolution = SplitConvolution<ScaleConvolution<FFTConvolution<T, FFTTag>>,FinegrainedPartitionnedFFTConvolution<T, FFTTag>>;


  template<typename A, typename B>
  struct PartitionAlgo< SplitConvolution<A,B> > {
    using NonAtomicConvolution = SplitConvolution<A,B>;
    using Delegate = typename NonAtomicConvolution::LateHandler;

    using SetupParam = typename Delegate::SetupParam;
    using PS = PartitionningSpec<SetupParam>;
    using PSpecs = PartitionningSpecs<SetupParam>;

    static PSpecs run(int n_channels, int n_audio_frames_per_cb, int size_impulse_response) {
      return PartitionAlgo<Delegate>::run(n_channels, n_audio_frames_per_cb, size_impulse_response);
    }
  };
}
