
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
   * Creates a 0-latency convolution by combining 'n' sub-convolutions
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
    }
    
    /*
     *   Unless the count of coefficients is of the form
     *     2^n - 1, some padding occurs.
     */
    // TODO we could require that 'A' tells what amount of underlying storage it will need, by number of coefficients,
    // and then this class could allocate the memory in one big chunk, and split it among 'A's.
    // This way, step / get could benefit from better memory locality, and memory accesses may be more predictable.
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
          // We pad up-to sizeBlock. Benchmarks showed that this is time-wise better
          // than padding to the next power of 2 + delaying input.
          auto nextPowOf2 = ceil_power_of_two(withPadding.size());
          withPadding.resize(sizeBlock);
          v[i].setCoefficients(withPadding);
        }
        else {
          v[i].setCoefficients({start,it});
        }
        if(v[i].getLatency() != sizeBlock - 1) {
          // This breaks the class logic, and would lead to wrong results.
          //   (for example, ScaleConvolution<FinegrainedPartitionnedFFTConvolution<float>> is not usable with this class.)
          LG(ERR,"ScaleConvolution is applied to a type that doesn't respect the latency constraints.");
          clear();
          return;
        }
      }
    }
    
    bool isValid() const {
      return std::all_of(v.begin(), v.end(), [](auto & e) -> bool { return e.isValid(); });
    }
    
    void step(FPT val) {
      std::for_each(v.begin(), v.end(), [val](auto & e) { e.step(val); });
    }
    
    auto get() const {
      FPT r{};
      for(auto const & e : v) {
        r += e.get();
      }
      return r;
      // TODO replace once supported by gcc
      //      return std::transform_reduce(v.begin(), v.end(), FPT{}, std::plus<>(), [](auto & e) { return e.get(); });
    }
    
    auto getEpsilon() const {
      FPT r{};
      for(auto const & e : v) {
        r += e.getEpsilon();
      }
      return r;
      // TODO replace once supported by gcc
      //      return std::transform_reduce(v.begin(), v.end(), FPT{}, std::plus<>(), [](auto & e) { return e.getEpsilon(); });
    }
    
    auto getLatency() const { return 0; }
    
  private:
    std::vector<A> v;
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
