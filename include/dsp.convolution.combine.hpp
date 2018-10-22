
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

  static constexpr int undefinedSplit = -1;
  static constexpr int noSplit = -2;

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

    static_assert(std::is_same_v<FPT,typename B::FPT>);

    bool isZero() const {
      return a.isZero() && b.isZero();
    }
    void reset() {
      a.reset();
      b.reset();
    }

    // The split should be equal to "latency of B - latency of A"
    // or negative to mean that no split exists.
    void setSplit(int s) {
      split = s;
    }

    void setCoefficients(a64::vector<FPT> coeffs_) {
      if(unlikely(undefinedSplit == split)) {
        throw std::logic_error("split was not set");
      }
      auto [rangeA,rangeB] = splitAt(split, coeffs_);
      a.setCoefficients(rangeA.materialize());
      if(rangeB.empty()) {
        b.reset();
        assert(b.isZero());
        return;
      }
      b.setCoefficients(rangeB.materialize());
      assert(!b.isZero());
      if(split + a.getLatency() == b.getLatency()) {
        return;
      }
      LG(ERR, "split : %d, a lat: %d, b lat: %d", split, a.getLatency(), b.getLatency());
      throw std::logic_error("inconsistent split (latencies are not adapted)");
    }

    bool isValid() const {
      return a.isValid() && (b.isValid() || b.isZero());
    }

    FPT step(FPT val) {
      auto res = a.step(val);
      if(!b.isZero()) {
        res += b.step(val);
      }
      return res;
    }
    
    auto debugStep(FPT val) {
      return std::make_pair(a.step(val), b.step(val));
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

  template<typename A, typename B>
  static std::ostream& operator<<(std::ostream& os, const typename SplitConvolution<A,B>::SetupParam& p)
  {
    using namespace std;
    os
    << "aParams:" << endl << p.aParams << endl
    << "bParams:" << endl << p.bParams << endl;
    return os;
  }

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

    bool isZero() const {
      return v.empty();
    }
    void reset() {
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
      reset();
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
      reset();
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
          reset();
          return;
        }
      }
      x.resize(pow2(n)); // including padding for biggest convolution
    }

    bool isValid() const {
      return std::all_of(v.begin(), v.end(), [](auto & e) -> bool { return e.isValid(); });
    }

    FPT step(FPT val) {
      auto it = v.begin();
      auto end = v.end();
      if(unlikely(it == end)) {
        return {};
      }

      x[progress] = typename RealSignal::value_type(val);
      ++progress;

      FPT r{};
      int nUpdates = 1 + (static_cast<int>(count_trailing_zeroes(progress))) - nDroppedConvolutions;
      if(nUpdates > 0) {
        typename RealSignal::const_iterator xBegin = x.begin();
        auto xBeginPadding = xBegin + progress;
        auto endUpdate = it + nUpdates;

        int lengthInputBeforePadding = pow2(nDroppedConvolutions);
        for(;it != endUpdate; ++it, lengthInputBeforePadding <<= 1) {
          // the start location should be 16 byte aligned for best performance.
          r += it->doStep(xBeginPadding - lengthInputBeforePadding);
          // and it's important that x[progress] to x[progress+lengthInputBeforePadding-1] are 0 (padding)
        }
      }

      for(; it!= end; ++it) {
        r += it->doStep();
      }

      if(nUpdates == v.size()) {
        Assert(2*progress == x.size());
        // fill with zeros so that padding is appropriate in the future.
        auto start = x.begin();
        auto end = start + progress; // not x.end() because by design, the second half of the vector is already zero-ed.
        std::fill(start, end, typename RealSignal::value_type(0));
        progress = 0;
      }

      return r;
    }
    double getEpsilon() const {
      return epsilonOfNaiveSummation(v);
    }

    auto getLatency() const {
      if(v.empty()) {
        return 0;
      }
      return v[0].getLatency();
    }

  private:
    std::vector<A> v;
    RealSignal x;
    unsigned int progress = 0;
    int nDroppedConvolutions = 0;
  };


  /*
   
   On filters, and how to chose them (here, filter and convolution means the same thing).

   Summary
   -------

   ********************************************************************************
   *
   * For ** live ** audio processing, use 'ZeroLatencyScaledFineGrainedPartitionnedConvolution':
   *   it has the smallest worst cost per audio callback.
   *
   * For ** offline ** audio processing, use 'OptimizedFIRFilter':
   *   it has the smallest average cost per sample.
   *
   * Both of them are 0-latency and are designed to scale both for very high and very low count of coefficients.
   ********************************************************************************

   Details
   -------

   There are 3 metrics that can be optimized:

   - average cost :
   How much CPU will I need, on average, to compute a single sample?
   This metric should be minimized when doing offline audio processing, so as to ensure that
   the task is executed as fast as possible.
   
   - "worst audiocallback" cost :
   How much CPU will I need, at worst, during a single audio callback, to compute the corresponding samples?
   Note that the size of the audio callback is important here.
   This metric should be optimized when doing live audio processing,
   to ensure that the audio callback meets its deadline.

   - latency :
   does my filter / convolution induce any latency?
   This metric should be optimized when doing live audio processing
   with an instrumentist playing a virtual instrument, because
   it's very hard / unpleasant to play an instrument that has noticeable latency.

   Here is the mapping from filter type to the kind of metric that is optimized by that filter:

   |----------------------------------------------------|----------------------------|-----------|
   |                                                    |     Optimal time cost      |           |
   |                                                    |----------------------------|           |
   |                                                    | on average | in worst case |           |
   | Filter type                                        | per sample | per callback  | 0-latency |
   |----------------------------------------------------|------------|---------------|-----------|
   | FIRFilter                                          |     .      |       .       |    X      |
   | OptimizedFIRFilter                                 |     X      |       .       |    X      |
   | FinegrainedPartitionnedFFTConvolution              |     .      |       X       |    .      |
   | ZeroLatencyScaledFineGrainedPartitionnedConvolution|     .      |       X       |    X      |
   |----------------------------------------------------|------------|---------------|-----------|

   */

  template<typename T, typename FFTTag = fft::Fastest>
  using OptimizedFIRFilter =
    SplitConvolution <
  // handle the first coefficients using brute force convolution (memory locality makes it faster than ScaleConvolution):
      FIRFilter<T>,
  // handle subsequent coefficients using FFTs, where each successive FFT doubles in size:
      ScaleConvolution <
        FFTConvolutionCore<T, FFTTag>
      >
    >;

  template<typename T, typename FFTTag = fft::Fastest>
  using ZeroLatencyScaledFineGrainedPartitionnedConvolution =
    SplitConvolution<
  // handle the first coefficients using a zero-latency filter:
      OptimizedFIRFilter<T, FFTTag>,
  // handle subsequent coefficients using a filter optimized for worst audio callback cost:
  SplitConvolution<
    FinegrainedPartitionnedFFTConvolution<T, FFTTag>, // full resolution
    SubSampled<
      LatencySemantic::DiracPeak,
      Delayed<

  SplitConvolution<
    FinegrainedPartitionnedFFTConvolution<T, FFTTag>, // half resolution
    SubSampled<
      LatencySemantic::DiracPeak,
      Delayed<

  SplitConvolution<
    FinegrainedPartitionnedFFTConvolution<T, FFTTag>, // quarter resolution
    SubSampled<
      LatencySemantic::DiracPeak,
      Delayed<

  FinegrainedPartitionnedFFTConvolution<T, FFTTag> // heighth resolution
  
      >
    >
  >
  
      >
    >
  >

      >
    >
  >

  >;


  /*
   Assuming that scales have resolutions of 1,2,4,8...
   for every scale, the number of blocks = { N(k), k in 1 .. S }
   N = sum (n in 0 .. (S-1)) 2^n N(k)
   and since we want the number of blocks to be almost equal, we will solve this:
   N = sum (n in 0 .. (S-1)) 2^n NF
   NF = N / (2^S - 1)
   */
  static inline int get_sz_for_single_scale(int response_sz, int n_scales) {
    int res = ceil(response_sz / static_cast<double>(pow2(n_scales) - 1));
    // must be even !
    if(res%2) {
      return res + 1;
    }
    return res;
  }

  template<typename T>
  struct EarlyestDeepest {
    using type = T;
  };

  template<typename A, typename B>
  struct EarlyestDeepest< SplitConvolution<A,B> > {
    using type = typename EarlyestDeepest<A>::type;
  };
  template<typename A>
  struct EarlyestDeepest< Delayed<A> > {
    using type = typename EarlyestDeepest<A>::type;
  };
  template<LatencySemantic Lat, typename A>
  struct EarlyestDeepest< SubSampled<Lat, A> > {
    using type = typename EarlyestDeepest<A>::type;
  };

  
  template<typename A, typename B>
  struct PartitionAlgo< SplitConvolution<A,B> > {
    using NonAtomicConvolution = SplitConvolution<A,B>;
    
  private:
    using LateHandler = typename EarlyestDeepest<typename NonAtomicConvolution::LateHandler>::type;
  public:
    using SetupParam = typename LateHandler::SetupParam;
    using PS = PartitionningSpec<SetupParam>;
    using PSpecs = PartitionningSpecs<SetupParam>;

    static PSpecs run(int n_channels, int n_audio_frames_per_cb, int total_response_size, int n_scales) {
      auto late_handler_response_size_for_partition_size =
      [total_response_size, n_scales](int partition_size) -> int {
        // we substract the latency (that many coefficients will be handled by the early coefficients handler).
        // it favors long partitions because we don't take into account
        // the cost of the early coefficient handler, but for long responses, where optimization matters,
        // it is negligible.
        auto n_coeffs_early_handler = LateHandler::getLatencyForPartitionSize(partition_size);
        auto late_response_sz = total_response_size - n_coeffs_early_handler;
        // return the size for a single scale.
        return get_sz_for_single_scale(late_response_sz, n_scales);
      };
      return PartitionAlgo<LateHandler>::run(n_channels,
                                             n_scales,
                                             n_audio_frames_per_cb,
                                             late_handler_response_size_for_partition_size);
    }
  };

  
  enum class AudioProcessing {
    Callback, // here, we want 0 latency and the smallest worst case cost per audio callback.
    Offline // here, we want the smallest averaged cost per sample.
  };
  
  namespace detail {
    template<typename T, AudioProcessing P>
    struct OptimalFilter_;
    
    template<typename T>
    struct OptimalFilter_<T, AudioProcessing::Callback> {
      using type = ZeroLatencyScaledFineGrainedPartitionnedConvolution<T>;
    };
    
    template<typename T>
    struct OptimalFilter_<T, AudioProcessing::Offline> {
      using type = OptimizedFIRFilter<T>;
    };
  }
  
  template<typename T, AudioProcessing P>
  using ZeroLatencyFilter = typename detail::OptimalFilter_<T,P>::type;
 

  template<typename Algo>
  void debugSteps(Algo & c, int nMax=20) {
    using namespace std;

    cout << "--" << endl;

    auto v = c.debugStep(1);
    int k = 0;
    cout << to_string(k) << " : " << to_string(v.first)<< " " << to_string(v.second) << endl;
    ++k;
    for(; k<nMax; ++k) {
      v = c.debugStep(0);
      cout << to_string(k) << " : " << to_string(v.first)<< " " << to_string(v.second) << endl;
    }
  }
}
