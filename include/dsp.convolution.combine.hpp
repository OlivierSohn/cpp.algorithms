
namespace imajuscule
{
  template<typename C>
  double epsilonOfNaiveSummation(C const & cont) {
    using FPT = typename std::remove_reference_t<decltype(cont[0])>::FPT;
    double err = {};
    // inner epsilons:
    for(auto const & c : cont) {
      err += c.getEpsilon();
    }
    // and there are cont.size() additions :
    err += cont.size() * std::numeric_limits<FPT>::epsilon();
    return err;
  }
  
  /*
   * Creates a convolution scheme by combining 2 convolution schemes,
   * where A handles ** early ** coefficients, and B handles ** late ** coefficients.
   */
  template<typename A, typename B>
  struct SplitConvolution {
    using FPT = typename A::FPT;
    using EarlyHandler = A;
    using LateHandler = B;
    
    struct SetupParam {
      using AParam = typename EarlyHandler::SetupParam;
      using BParam = typename LateHandler::SetupParam;
      AParam aParams;
      BParam bParams;
    };
    
    static constexpr int undefinedSplit = -1;

    static_assert(std::is_same_v<FPT,typename B::FPT>);

    bool empty() const {
      return a.empty() && b.empty();
    }
    void clear() {
      a.clear();
      b.clear();
    }

    void setSplit(int s) {
      if(unlikely(split != undefinedSplit)) {
        throw std::logic_error("split set twice");
      }
      split = s;
    }

    void setCoefficients(a64::vector<FPT> coeffs_) {
      if(unlikely(undefinedSplit == split)) {
        throw std::logic_error("split was not set");
      }
      auto [rangeA,rangeB] = splitAt(split, coeffs_);
      a.setCoefficients(rangeA.materialize());
      b.setCoefficients(rangeB.materialize());
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

    double getEpsilon() const {
      return a.getEpsilon() + b.getEpsilon() + std::numeric_limits<FPT>::epsilon();
    }

    auto getLatency() const {
      return a.getLatency();
    }
    
    auto & getA() { return a; }
    auto & getB() { return b; }

  private:
    int split = undefinedSplit; // we have 'split' ** early ** coefficients, the rest are ** late ** coefficients.
    A a;
    B b;
  };

  /*
   * Creates a 0-latency convolution by combining 'n' sub-convolutions
   * of latencies 2^k-1, where k is in [0,n-1].
   */
  // The code is commented out because it is superseeded by 'ScaleConvolution' which has more efficient memory usage.
  // The code was not deleted because the simplicity of the implementation could inspire other implementations.
  /*
  template<typename A>
  struct ScaleConvolutionLegacy {
    using FPT = typename A::FPT;
    
    bool empty() const {
      return v.empty() || std::all_of(v.begin(), v.end(), [](auto & e) -> bool { return e.empty(); });
    }
    void clear() {
      v.clear();
    }
    
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
          withPadding.resize(sizeBlock);
          v[i].setCoefficients(withPadding);
        }
        else {
          v[i].setCoefficients({start,it});
        }
        if(v[i].getLatency() != sizeBlock - 1) {
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
    }
    double getEpsilon() const {
      return epsilonOfNaiveSummation(v);
    }
    
    auto getLatency() const { return 0; }
    
  private:
    std::vector<A> v;
  };
*/


  /*
   * Creates a 0-latency convolution by combining 'n' sub-convolutions
   * of latencies 2^k-1, where k is in [0,n-1].
   *
   * The first convolutions can be dropped (see the constructor).
   */
  template<typename A>
  struct ScaleConvolution {
    using FPT = typename A::FPT;
    using RealSignal = typename A::RealSignal;
    
    struct SetupParam {
      int nDropped = 5; // see 'TEST(BenchmarkConvolutions, scaled_vs_brute)' :
      // 5 gives the fastest results on OSX / intel core i7.
    };
    
    bool empty() const {
      return v.empty() || std::all_of(v.begin(), v.end(), [](auto & e) -> bool { return e.empty(); });
    }
    void clear() {
      v.clear();
      x.clear();
      progress = 0;
    }
    
    /*
     * Sets the number of early convolutions that are dropped.
     * The coefficients associated to these dropped convolutions are typically handled
     * by another convolution handler. (there are '2^nDropped - 1' such coefficients)
     *
     * Also calls 'clear'. Hence, if you want to set the number of dropped convolutions,
     * you should do so before calling 'setCoefficients'.
     *
     * Note that the default value is optimal for the system I developped on.
     * see 'TEST(BenchmarkConvolutions, scaled_vs_brute)' on how to determine the best value for a goven system.
     *
     * TODO ideally this value should be a global, computed at initialization time.
     */
    void applySetup(SetupParam const & p) {
      clear();
      nDroppedConvolutions = p.nDropped;
    }
    
    int getEarlyDroppedConvolutions() const {
      return nDroppedConvolutions;
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
      // note that these dropped coefficients are not passed to this function
      auto nFirstCoefficients = pow2(nDroppedConvolutions)-1;
      auto n = power_of_two_exponent(ceil_power_of_two(nFirstCoefficients+s+1));
      Assert(nDroppedConvolutions < n);
      v.resize(n-nDroppedConvolutions);
      auto it = coeffs_.begin();
      auto end = coeffs_.end();
      for(int i=nDroppedConvolutions; i<n; ++i) {
        auto & conv = v[i-nDroppedConvolutions];
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
          conv.setCoefficients2(withPadding);
        }
        else {
          conv.setCoefficients2({start,it});
        }
        if(conv.getLatency() != sizeBlock - 1) {
          // This breaks the class logic, and would lead to wrong results.
          //   (for example, ScaleConvolution<FinegrainedPartitionnedFFTConvolution<float>> is not usable with this class.)
          LG(ERR,"ScaleConvolution is applied to a type that doesn't respect the latency constraints.");
          clear();
          return;
        }
      }
      x.resize(pow2(n)); // including padding for biggest convolution
    }
    
    bool isValid() const {
      return std::all_of(v.begin(), v.end(), [](auto & e) -> bool { return e.isValid(); });
    }
    
    void step(FPT val) {
      auto it = v.begin();
      auto end = v.end();
      if(unlikely(it == end)) {
        return;
      }

      x[progress] = typename RealSignal::value_type(val);
      ++progress;

      int nUpdates = 1 + (static_cast<int>(count_trailing_zeroes(progress))) - nDroppedConvolutions;
      if(nUpdates > 0) {
        typename RealSignal::const_iterator xBegin = x.begin();
        auto xBeginPadding = xBegin + progress;
        auto endUpdate = it + nUpdates;
        
        int lengthInputBeforePadding = pow2(nDroppedConvolutions);
        for(;it != endUpdate; ++it, lengthInputBeforePadding <<= 1) {
          // the start location should be 16 byte aligned for best performance.
          it->doStep(xBeginPadding - lengthInputBeforePadding);
          // and it's important that x[progress] to x[progress+lengthInputBeforePadding-1] are 0 (padding)
        }
      }

      for(; it!= end; ++it) {
        it->doStep();
      }
      
      if(nUpdates == v.size()) {
        Assert(2*progress == x.size());
        // fill with zeros so that padding is appropriate in the future.
        auto start = x.begin();
        auto end = start + progress; // not x.end() because by design, the second half of the vector is already zero-ed.
        std::fill(start, end, typename RealSignal::value_type(0));
        progress = 0;
      }
    }
    
    auto get() const {
      FPT r{};
      for(auto const & e : v) {
        r += e.get();
      }
      return r;
    }
    double getEpsilon() const {
      return epsilonOfNaiveSummation(v);
    }
    
    auto getLatency() const { return 0; }
    
  private:
    std::vector<A> v;
    RealSignal x;
    unsigned int progress = 0;
    int nDroppedConvolutions = 0;
  };
 
  template<typename C>
  void setPartitionSize(C & c, int sz) {
    c.getB().set_partition_size(sz);
    c.setSplit( c.getB().getLatency() - c.getA().getLatency() );
  }

  
  template<typename C, typename SetupParam>
  void applySetup(C&c, SetupParam const & p) {
    c.applySetup(p);
  }
  
  template<typename A, typename B, typename SetupParam>
  void applySetup(SplitConvolution<A,B> &c, SetupParam const & p) {
    applySetup(c.getA(), p.aParams);
    applySetup(c.getB(), p.bParams);
  }
  template<typename A, typename B, typename SetupParam>
  void applySetup(SplitConvolution<A,ScaleConvolution<B>> &c, SetupParam const & p) {
    applySetup(c.getA(), p.aParams);
    applySetup(c.getB(), p.bParams);
    c.setSplit(pow2(c.getB().getEarlyDroppedConvolutions()) - 1);
  }


  template<typename T>
  using ZeroLatencyBruteConvolution = FIRFilter<T>;

  template<typename T, typename FFTTag = fft::Fastest>
  using ZeroLatencyBruteFineGrainedPartitionnedConvolution = SplitConvolution<ZeroLatencyBruteConvolution<T>,FinegrainedPartitionnedFFTConvolution<T, FFTTag>>;

  template<typename T, typename FFTTag = fft::Fastest>
  using ZeroLatencyScaledFineGrainedPartitionnedConvolution =
    SplitConvolution<
      SplitConvolution<
        FIRFilter<T>,
        ScaleConvolution<
          FFTConvolutionCore<T, FFTTag>
        >
      >,
      FinegrainedPartitionnedFFTConvolution<T, FFTTag>
  >;


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
