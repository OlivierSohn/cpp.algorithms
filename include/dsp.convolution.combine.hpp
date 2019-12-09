
namespace imajuscule {

template<std::size_t i, std::size_t j>
auto arrayConcat(std::array<int, i> const & a,
                 std::array<int, j> const & b) -> std::array<int, i+j>
{
    std::array<int, i+j> r;
    for(int k =0; k<i; ++k) {
        r[k] = a[k];
    }
    for(int k =0; k<j; ++k) {
        r[i+k] = b[k];
    }
    return r;
}

template<std::size_t i, std::size_t j>
auto arraySplit(std::array<int, i+j> const & r) -> std::pair<std::array<int, i>, std::array<int, j>>
{
    std::array<int, i> a;
    std::array<int, j> b;
    for(int k =0; k<i; ++k) {
        a[k] = r[k];
    }
    for(int k =0; k<j; ++k) {
        b[k] = r[i+k];
    }
    return {a,b};
}

template<typename T>
struct pointed {
    using type = T;
    static auto const & get(T const & t) {
        return t;
    }
};
template<typename T>
struct pointed<std::unique_ptr<T>> {
    using type = T;
    static auto const & get(std::unique_ptr<T> const & t) {
        return *(t.get());
    }
};

template<typename C>
double epsilonOfNaiveSummation(C const & cont) {
    using Pointed = pointed<typename C::value_type>;
    using FPT = typename std::remove_reference_t<typename Pointed::type>::FPT;
    double err = {};
    // inner epsilons:
    for(auto const & c : cont) {
      err += Pointed::get(c).getEpsilon();
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
    static constexpr int nComputePhaseable = A::nComputePhaseable + B::nComputePhaseable;
    static constexpr int nCoefficientsFadeIn = A::nCoefficientsFadeIn;
    static constexpr bool has_subsampling = LateHandler::has_subsampling; // we 'could' take earlyhandler into account, too.

    struct SetupParam : public Cost {
      using AParam = typename EarlyHandler::SetupParam;
      using BParam = typename LateHandler::SetupParam;
      AParam aParams;
      BParam bParams;
        
        SetupParam(AParam a, BParam b)
        : aParams(a)
        , bParams(b)
        {}

        void logSubReport(std::ostream & os) const override {
            aParams.logSubReport(os);
            bParams.logSubReport(os);
        }
    };
      
      void logComputeState(std::ostream & os) const {
          os << "[0.." << split-1 << "]" << std::endl;
          {
              IndentingOStreambuf indent{std::cout};
              a.logComputeState(os);
          }
          os << "[" << split << "..]" << std::endl;
          {
              IndentingOStreambuf indent{std::cout};
              b.logComputeState(os);
          }
      }

    static_assert(std::is_same_v<FPT,typename B::FPT>);

    bool isZero() const {
      return a.isZero() && b.isZero();
    }
    void reset() {
      a.reset();
      b.reset();
      split = undefinedSplit;
    }
    void flushToSilence() {
      a.flushToSilence();
      b.flushToSilence();
    }
            
      void setup(const SetupParam & p) {
          a.setup(p.aParams);
          b.setup(p.bParams);

          if(!b.isValid()) { // for case where lower resolution tail is not used
            split = noSplit;
          }
          else {
            split = B::nCoefficientsFadeIn + b.getLatency() - a.getLatency();
          }
      }

    void setCoefficients(a64::vector<FPT> coeffs_) {
      if(unlikely(undefinedSplit == split)) {
        throw std::logic_error("undefined split");
      }
      auto [rangeA,rangeB] = splitAt(split, coeffs_);
      
      if(rangeB.empty()) {
        a.setCoefficients(rangeA.materialize());
        b.reset();
        assert(b.isZero());
        return;
      }
      else {
        // prepend overlapping coefficients to b range:
        rangeB.setBegin(rangeB.begin() - B::nCoefficientsFadeIn);
        if(rangeB.begin() < rangeA.begin()) {
          LG(ERR, "sz_coeffs: %d, split: %d, nFadeIn : %d", coeffs_.size(), split, B::nCoefficientsFadeIn);
          throw std::logic_error("cannot cross-fade the coefficients");
        }
        a.setCoefficients(withLinearFadeOut(B::nCoefficientsFadeIn,rangeA.materialize()));
        b.setCoefficients(withLinearFadeIn( B::nCoefficientsFadeIn,rangeB.materialize()));
      }
      assert(!b.isZero());
      if(split + a.getLatency() == B::nCoefficientsFadeIn + b.getLatency()) {
        return;
      }
      LG(ERR, "split : %d, a lat: %d, b lat: %d, b init: %d", split, a.getLatency(), b.getLatency(), B::nCoefficientsFadeIn);
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
      template<typename FPT2>
      void stepAddVectorized(FPT2 const * const input_buffer,
                             FPT2 * output_buffer,
                             int nSamples)
      {
          a.stepAddVectorized(input_buffer,
                              output_buffer,
                              nSamples);
          if(!b.isZero()) {
              b.stepAddVectorized(input_buffer,
                                  output_buffer,
                                  nSamples);
          }
      }
      template<typename FPT2>
      void stepAssignVectorized(FPT2 const * const input_buffer,
                                FPT2 * output_buffer,
                                int nSamples)
      {
          a.stepAssignVectorized(input_buffer,
                                 output_buffer,
                                 nSamples);
          if(!b.isZero()) {
              b.stepAddVectorized(input_buffer,
                                  output_buffer,
                                  nSamples);
          }
      }
      
      template<typename FPT2>
      void stepAddInputZeroVectorized(FPT2 * output_buffer,
                                      int nSamples)
      {
          a.stepAddInputZeroVectorized(output_buffer,
                                       nSamples);
          if(!b.isZero()) {
              b.stepAddInputZeroVectorized(output_buffer,
                                           nSamples);
          }
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

      std::array<int, nComputePhaseable> getComputePeriodicities() const {
          return arrayConcat(a.getComputePeriodicities(),
                             b.getComputePeriodicities());
      }
      // in [0, getComputePeriodicity())
      std::array<int, nComputePhaseable> getComputeProgresses() const {
        return arrayConcat(a.getComputeProgresses(),
                           b.getComputeProgresses());
      }
      void setComputeProgresses(std::array<int, nComputePhaseable> const & progresses) {
          auto [pa,pb] = arraySplit<A::nComputePhaseable,B::nComputePhaseable>(progresses);
          a.setComputeProgresses(pa);
          if(!b.isZero()) {
              b.setComputeProgresses(pb);
          }
      }

    auto & getA() { return a; }
    auto & getB() { return b; }
      auto const & getA() const { return a; }
      auto const & getB() const { return b; }

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

      
      struct CountDroppedScales {
          constexpr explicit CountDroppedScales(int n)
          : n(n)
          {}

          constexpr int toInteger() const {
            return n;
          }
      private:
            int n;
      };

struct ScaleConvolution_ {
    static constexpr int latencyForDroppedConvolutions(CountDroppedScales nDropped) {
        return static_cast<int>(pow2(static_cast<size_t>(nDropped.toInteger())))-1;
    }

    /*
     * The default value is optimal for the system I developped on (OSX / Macbook Air 2015 / intel core i7).
     * see 'TEST(BenchmarkConvolutions, scaled_vs_brute)' on how to determine the best value for a given system.
     *
     * TODO ideally this value should be a global, computed at initialization time.
     */
    static constexpr CountDroppedScales nDroppedOptimalFor_Split_Bruteforce_Fft = CountDroppedScales(6);
};

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
    using Tag = typename A::Tag;
    static constexpr auto zero_signal = fft::RealSignal_<Tag, FPT>::zero;

    static constexpr int nComputePhaseable = 1;
    static constexpr int nCoefficientsFadeIn = 0;
    static_assert(A::nCoefficientsFadeIn == 0); // else we need to handle them

    // note that it wouldn't make much sense for 'A' to have subsampling in a ScaleConvolution
      static constexpr bool has_subsampling = A::has_subsampling;
      
    struct SetupParam : public Cost {
        SetupParam(CountDroppedScales nDropped = ScaleConvolution_::nDroppedOptimalFor_Split_Bruteforce_Fft)
        : nDropped(nDropped)
        {}
        
      /*
       * The number of early convolutions that are dropped.
       * The coefficients associated to these dropped convolutions are typically handled
       * by another convolution handler. (there are '2^nDropped - 1' such coefficients)
       *
       * (todo refactor to remove this warning : every PartitionAlgo should explicitely set this value)
       * WARNING: Do not change the default value unless you know what you are doing:
       * PartitionAlgo< ZeroLatencyScaledFineGrainedPartitionnedConvolution<T,FFTTag> >
       * assumes that this value will be the default value.
       */
      CountDroppedScales nDropped;
      
      int getImpliedLatency() const {
        return ScaleConvolution_::latencyForDroppedConvolutions(nDropped);
      }
      
      void logSubReport(std::ostream & os) const override {
        os << "Scaling, dropped: " << nDropped.toInteger() << std::endl;
      }
    };
      
      void logComputeState(std::ostream & os) const {
          os << "Scaling ["<< progress <<"/"<< x.size()/2 <<"], dropped: " << nDroppedConvolutions.toInteger() << std::endl;
          IndentingOStreambuf indent(os);
          int i=nDroppedConvolutions.toInteger();
          for(auto const & algo : v)
          {
              ++i;
              os << i << std::endl;
              IndentingOStreambuf indent2(os);
              algo.logComputeState(os);
          }
      }
      
    bool isZero() const {
      return v.empty();
    }
    void reset() {
      v.clear();
      x.clear();
      progress = 0;
      nDroppedConvolutions = CountDroppedScales(0);
    }
    void flushToSilence() {
      for(auto & c:v) {
        c.flushToSilence();
      }
      if(!x.empty()) {
        zero_signal(x);
      }
      progress = 0;
    }

    void setup(SetupParam const & p) {
      reset();
      nDroppedConvolutions = p.nDropped;
    }
      
    /*
     *   Unless the count of coefficients is of the form
     *     2^n - 1, some padding occurs.
     */
    // TODO we could require that 'A' tells what amount of underlying storage it will need, by number of coefficients,
    // and then this class could allocate the memory in one big chunk, and split it among 'A's.
    // This way, step / get could benefit from better memory locality, and memory accesses may be more predictable.
    void setCoefficients(a64::vector<FPT> coeffs_) {
      auto s = coeffs_.size();
      assert( s > 0 );
      // note that these dropped coefficients are not passed to this function
      auto nFirstCoefficients = getLatency();
      auto n = power_of_two_exponent(ceil_power_of_two(nFirstCoefficients+s+1));
      Assert(nDroppedConvolutions.toInteger() < n);
      v.resize(n-nDroppedConvolutions.toInteger());
      auto it = coeffs_.begin();
      auto end = coeffs_.end();
      for(int i=nDroppedConvolutions.toInteger(); i<n; ++i) {
        auto & conv = v[i-nDroppedConvolutions.toInteger()];
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
          v.clear();
          return;
        }
      }
      x.resize(pow2(n)); // including padding for biggest convolution
    }

    bool isValid() const {
        if(nDroppedConvolutions.toInteger() < 0) {
            return false;
        }
        if(v.empty()) {
            return true;
        }
      return std::all_of(v.begin(), v.end(), [](auto & e) -> bool { return e.isValid(); });
    }

    FPT step(FPT val) {
      if(unlikely(isZero())) {
        return {};
      }
      return doStep(val);
    }
  private:
    FPT doStep(FPT val) {
      auto it = v.begin();
      auto end = v.end();

      x[progress] = typename RealSignal::value_type(val);
      ++progress;

      FPT r{};
      int nUpdates = 1 + (static_cast<int>(count_trailing_zeroes(progress))) - nDroppedConvolutions.toInteger();
      if(nUpdates > 0) {
        typename RealSignal::const_iterator xBegin = x.begin();
        auto xBeginPadding = xBegin + progress;
        auto endUpdate = it + nUpdates;

        int lengthInputBeforePadding = pow2(nDroppedConvolutions.toInteger());
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
public:
      template<typename FPT2>
      void stepAddVectorized(FPT2 const * const input_buffer,
                             FPT2 * output_buffer,
                             int nSamples)
      {
          if(unlikely(isZero())) {
            return;
          }
          for(int i=0; i<nSamples; ++i) {
              output_buffer[i] += doStep(input_buffer[i]);
          }
      }
      template<typename FPT2>
      void stepAssignVectorized(FPT2 const * const input_buffer,
                                FPT2 * output_buffer,
                                int nSamples)
      {
          if(unlikely(isZero())) {
            return;
          }
          for(int i=0; i<nSamples; ++i) {
              output_buffer[i] = doStep(input_buffer[i]);
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
              output_buffer[i] += doStep({});
          }
      }
      
    double getEpsilon() const {
      return epsilonOfNaiveSummation(v);
    }

    auto getLatency() const {
        return ScaleConvolution_::latencyForDroppedConvolutions(nDroppedConvolutions);
    }
      std::array<int, 1> getComputePeriodicities() const {
          int res = x.size() / 2;
          Assert(res == getBiggestScale());
        return {res};
      }
      // in [0, getComputePeriodicity())
      std::array<int, 1> getComputeProgresses() const {
        return {static_cast<int>(progress)};
      }
      void setComputeProgresses(std::array<int, 1> const & progresses) {
        auto const p = progresses[0];
        while(getComputeProgresses()[0] != p) {
            step(0);
        }
      }
      
      int countCoefficients() const {
          return pow2(nDroppedConvolutions + v.size()) - pow2(nDroppedConvolutions);
      }
      
      int getBiggestScale() const {
          if(v.empty()) {
              return 0;
          }
          return pow2(nDroppedConvolutions.toInteger() + v.size() - 1);
      }

  private:
    std::vector<A> v;
    RealSignal x;
    unsigned int progress = 0;
    CountDroppedScales nDroppedConvolutions = CountDroppedScales(0);
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
  
  template<typename T, typename FFTTag = fft::Fastest, LatencySemantic Lat = LatencySemantic::DiracPeak>
  using ZeroLatencyScaledFineGrainedPartitionnedSubsampledConvolution =
    SplitConvolution<
  // handle the first coefficients using a zero-latency filter:
      OptimizedFIRFilter<T, FFTTag>,
  // handle subsequent coefficients using a filter optimized for worst audio callback cost:
  SplitConvolution<
    FinegrainedPartitionnedFFTConvolution<T, FFTTag>, // full resolution
    SubSampled<
      Lat,
      Delayed<

  SplitConvolution<
    FinegrainedPartitionnedFFTConvolution<T, FFTTag>, // half resolution
    SubSampled<
      Lat,
      Delayed<

  SplitConvolution<
    FinegrainedPartitionnedFFTConvolution<T, FFTTag>, // quarter resolution
    SubSampled<
      Lat,
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

      template<typename T, typename FFTTag = fft::Fastest>
      using ZeroLatencyScaledFineGrainedPartitionnedConvolution =
        SplitConvolution<
      // handle the first coefficients using a zero-latency filter:
          OptimizedFIRFilter<T, FFTTag>,
      // handle subsequent coefficients using a filter optimized for worst audio callback cost:
          FinegrainedPartitionnedFFTConvolution<T, FFTTag>
      >;
      template<typename T, typename FFTTag = fft::Fastest>
      using ZeroLatencyScaledAsyncConvolution =
      SplitConvolution<
      // handle the first coefficients using a zero-latency filter:
        OptimizedFIRFilter<T, FFTTag>,
      // handle subsequent coefficients using an asynchronous scaling convolution
        AsyncCPUConvolution < ScaleConvolution < FFTConvolutionCore<T, FFTTag> > >
      >;      
      template<typename T, typename FFTTag = fft::Fastest>
      using ZeroLatencyScaledAsyncConvolutionOptimized =
      SplitConvolution<
      // handle the first coefficients using
        ZeroLatencyScaledFineGrainedPartitionnedConvolution<T, FFTTag>,
      // handle subsequent coefficients using an asynchronous scaling convolution
        AsyncCPUConvolution < ScaleConvolution < FFTConvolutionCore<T, FFTTag> > >
      >;

      template<typename T, typename FFTTag>
      struct PartitionAlgo< OptimizedFIRFilter<T, FFTTag> > {
          using Convolution = OptimizedFIRFilter<T, FFTTag>;
          using SetupParam = typename Convolution::SetupParam;
          using PS = PartitionningSpec<SetupParam>;
          using PSpecs = PartitionningSpecs<SetupParam>;
          
          static PSpecs run(int n_channels,
                            int n_audio_channels,
                            int n_audio_frames_per_cb,
                            int total_response_size,
                            double frame_rate,
                            double max_avg_time_per_sample,
                            std::ostream & os) {
              // there is no variable to optimize with OptimizedFIRFilter:
              PS minimalPs;
              minimalPs.cost = SetupParam({},{});
              minimalPs.cost->setCost(0.);
              return {
                  minimalPs,
                  minimalPs
              };
          }
      };

      enum class SimulationPhasingMode {
        On,
        Off
      };
      struct SimulationPhasing {
        SimulationPhasingMode mode;
        std::optional<int> groupSize;
      
          static SimulationPhasing no_phasing() {
            return {SimulationPhasingMode::Off, {}};
          }
          static SimulationPhasing phasing_with_group_size(int sz) {
            return {SimulationPhasingMode::On, {sz}};
          }
      };
      
      static inline std::ostream & operator << (std::ostream & o, SimulationPhasing const & p) {
          if(p.mode == SimulationPhasingMode::Off) {
            o << "off";
          }
          else {
            o << "with group size : " << *p.groupSize;
          }
        return o;
      }


        template<typename Algo>
        struct PartitionAlgo< AsyncCPUConvolution<Algo> > {
            using Convolution = AsyncCPUConvolution<Algo>;
            using SetupParam = typename Convolution::SetupParam;
            using PS = PartitionningSpec<SetupParam>;
            using PSpecs = PartitionningSpecs<SetupParam>;

            struct Metrics {
                range<int> resultQueueEltsCountRange;
                
                // We are interested in the minimum size, by period, of the result queue.
                // If the minimum size keeps being smaller over the periods, it meens that
                // the background thread is too slow.

                double minResultQueueEltsCountRange_gradient = 0; // diff between last and first period over number of periods
                // a strictly negative gradient means that the system is not fast enough.
                
                int minResultQueueEltsCountRange_min_local_gradient = 0; // diff between consecutive periods
                
                void mergeWorst(Metrics const & o) {
                    resultQueueEltsCountRange.extend(o.resultQueueEltsCountRange);
                    minResultQueueEltsCountRange_gradient = std::min(minResultQueueEltsCountRange_gradient,
                                                     o.minResultQueueEltsCountRange_gradient);
                    minResultQueueEltsCountRange_min_local_gradient = std::min(minResultQueueEltsCountRange_min_local_gradient,
                                                               o.minResultQueueEltsCountRange_min_local_gradient);
                }
            };
            
            struct QueueMetrics {
                QueueMetrics(int period, int nPeriods)
                : nPeriods(nPeriods)
                , period_frames(period)
                {
                    if(nPeriods < 2) {
                        throw std::logic_error("2 periods are required to compute a gradient");
                    }
                    resultQueueEltsCountRangeByPeriod.reserve(nPeriods);
                    resultQueueEltsCountRangeByPeriod.resize(1);
                    signalQueueEltsCountRangeByPeriod.reserve(nPeriods);
                    signalQueueEltsCountRangeByPeriod.resize(1);
                }
                
            private:
                int const nPeriods;
                int const period_frames;
                int period = 0;
                int n_cur_frame = 0;
                bool hasErrorWorkerTooSlow = false;
                
                std::vector<range<int>> resultQueueEltsCountRangeByPeriod, signalQueueEltsCountRangeByPeriod;
            public:
                bool recordQueueSize(int resultQueueEltsCount,
                                     int signalQueueEltsCount,
                                     int countErrorsWorkerTooSlow,
                                     int nFrames) {
                    assert(period < nPeriods);
                    if(countErrorsWorkerTooSlow) {
                        hasErrorWorkerTooSlow = true;
                        return false;
                    }
                    n_cur_frame += nFrames;
                    if(n_cur_frame >= period_frames) {
                        n_cur_frame -= period_frames;
                        ++period;
                        if(period >= nPeriods) {
                            return false;
                        }
                        Assert(resultQueueEltsCountRangeByPeriod.size() < resultQueueEltsCountRangeByPeriod.capacity());
                        Assert(signalQueueEltsCountRangeByPeriod.size() < signalQueueEltsCountRangeByPeriod.capacity());
                        resultQueueEltsCountRangeByPeriod.emplace_back();
                        signalQueueEltsCountRangeByPeriod.emplace_back();
                    }
                    resultQueueEltsCountRangeByPeriod.back().extend(resultQueueEltsCount);
                    signalQueueEltsCountRangeByPeriod.back().extend(signalQueueEltsCount);
                    return true;
                }
                
                std::vector<range<int>> const & getResultQueueEltsCountRangeByPeriod () const {
                    return resultQueueEltsCountRangeByPeriod;
                }
                std::vector<range<int>> const & getSignalQueueEltsCountRangeByPeriod () const {
                    return signalQueueEltsCountRangeByPeriod;
                }

                std::optional<Metrics> getMetrics() const {
                    if(hasErrorWorkerTooSlow) {
                      return {};
                    }
                    Metrics m;
                    for(auto const & r : resultQueueEltsCountRangeByPeriod) {
                        m.resultQueueEltsCountRange.extend(r);
                    }
                    {
                        std::optional<int> prevMinQueueSize;
                        std::optional<int> minGradient;
                        bool first = true;
                        for(auto const & r : resultQueueEltsCountRangeByPeriod) {
                            if(first) {
                                // discard first period for gradient computation
                                // (we are not in a stable mode at the beginning of the first period)
                                first = false;
                                break;
                            }
                            auto minResultQueueSizeDuringPeriod = r.getMin();
                            if(!prevMinQueueSize) {
                                prevMinQueueSize = minResultQueueSizeDuringPeriod;
                            }
                            else {
                                auto gradient = minResultQueueSizeDuringPeriod - *prevMinQueueSize;
                                if(!minGradient || *minGradient > gradient) {
                                    minGradient = gradient;
                                }
                            }
                        }
                        m.minResultQueueEltsCountRange_min_local_gradient = *minGradient;
                    }
                    if(resultQueueEltsCountRangeByPeriod.size() >= 3)
                    {
                        auto & first = *(resultQueueEltsCountRangeByPeriod.begin() + 1);
                        auto & last = resultQueueEltsCountRangeByPeriod.back();
                        m.minResultQueueEltsCountRange_gradient = (last.getMin() - first.getMin()) / static_cast<double>(nPeriods);
                    }
                    else {
                        throw std::logic_error("not enough periods to compute a gradient");
                    }
                    return m;
                }
            };

            static constexpr int numPeriods = 5; // the first period is not taken into account for gradient computation

            // todo use a better approach than this threshold
            static constexpr double normalizedGradientThreshold = -0.2;

            static PSpecs run(int n_channels,
                              int n_audio_channels,
                              int n_audio_frames_per_cb,
                              int total_response_size,
                              double frame_rate,
                              double max_avg_time_per_sample,
                              std::ostream & os,
                              SimulationPhasing const & phasing,
                              CountDroppedScales const & nAsyncScalesDropped) // todo gerer en args...
            {
                auto res = computeQueueMetrics(n_channels,
                                                n_audio_channels,
                                                n_audio_frames_per_cb,
                                                total_response_size,
                                                frame_rate,
                                                nAsyncScalesDropped,
                                                phasing);
                if(!res) {
                    os << "Synchronous thread is too slow" << std::endl;
                    return {};
                }
                auto [periodically, queueMetrics] = *res;

                // If there was jitter during the simulation, the metrics might be pessimistic.
                // The jitter of the simulation could be analyzed using 'periodically'
                int const submission_period = n_audio_frames_per_cb;
                int minNumFramesPerQueue = 0;
                if(queueMetrics.empty()) {
                  // no late coeffs
                }
                else {
                    std::vector<Metrics> metrics;
                    metrics.reserve(queueMetrics.size());
                    for(auto const & m : queueMetrics) {
                        auto mayMetric = m.getMetrics();
                        if(!mayMetric) {
                            os << "Background worker is too slow : a result queue was empty" << std::endl;
                            return {};
                        }
                        metrics.push_back(*mayMetric);
                    }
                    
                    Metrics worst;
                    for(auto const & m : metrics) {
                      worst.mergeWorst(m);
                    }
                    
                    auto normalizedGradient = worst.minResultQueueEltsCountRange_gradient / (1 + worst.resultQueueEltsCountRange.getSpan());
                    if(normalizedGradient < normalizedGradientThreshold) {
                        os << "Background worker is too slow : normalizedGradient on min queue elts count = " << normalizedGradient << std::endl;
                        int i=0;
                        for(auto const & m : queueMetrics) {
                            std::cout << i << std::endl;
                            for(auto const & r : m.getResultQueueEltsCountRangeByPeriod()) {
                              std::cout << r.getMin() << " " << r.getMax() << std::endl;
                            }
                            std::cout << std::endl;
                            for(auto const & r : m.getSignalQueueEltsCountRangeByPeriod()) {
                              std::cout << r.getMin() << " " << r.getMax() << std::endl;
                            }
                            std::cout << std::endl;
                            ++i;
                        }
                        return {};
                    }
                    
                    constexpr int safeFactor = 4;
                    
                    minNumFramesPerQueue =
                    safeFactor * // to account for a system that is under heavier load
                    n_audio_frames_per_cb * // assumes that the submission_period was n_audio_frames_per_cb during the simulation
                    (1 +
                     worst.resultQueueEltsCountRange.getSpan());
                    // to avoid audio dropouts in case of the occurence of the race condition
                    // commented in the code of the worker of AsyncCPUConvolution:
                    minNumFramesPerQueue += n_audio_frames_per_cb;
                }
                os << "minNumFramesPerQueue=" << minNumFramesPerQueue << std::endl;

                int const queueSize = minNumFramesPerQueue/submission_period;
                
                PS minimalPs;
                minimalPs.cost = SetupParam{
                    submission_period,
                    queueSize,
                    {
                        nAsyncScalesDropped // todo il faudrait peut-etre revoir cela a la hausse, (quitte a augmenter la latence du late) poru reduire la somme total de calculs
                    }
                };
                minimalPs.cost->setCost(0.);
                return {
                    minimalPs,
                    minimalPs
                };
            }
        private:
            static std::optional<std::pair<std::optional<Periodically>,std::vector<QueueMetrics>>>
            computeQueueMetrics(int n_channels,
                                int n_audio_channels,
                                int n_audio_frames_per_cb,
                                int total_response_size,
                                double frame_rate,
                                CountDroppedScales const & nAsyncScalesDropped,
                                SimulationPhasing const & phasing)
            {
                std::vector<std::unique_ptr<Convolution>> async_convs;
                async_convs.reserve(n_channels);
                for(int i=0; i<n_channels; ++i) {
                    async_convs.push_back(std::make_unique<Convolution>());
                }

                // we use a huge queue so that we can have room for queue size variations
                int const num_frames_in_queue = std::max(10000, total_response_size/2);
        
                int const submission_period = n_audio_frames_per_cb;
                int queue_size = num_frames_in_queue / submission_period;
                if(queue_size * submission_period < num_frames_in_queue) {
                    ++queue_size;
                }
                
                // this is pessimistic, it will be more depending on the computed queue size
                // (the bigger the queue size,
                //  the bigger the latency of the latehandler,
                //  the more coefficients are handled by the early handler)
                int const nMinDroppedCoeffs = n_audio_frames_per_cb;
                int const nLateCoeffs = total_response_size-nMinDroppedCoeffs;
                if(nLateCoeffs <= 0) {
                    // the size of the queue doesn't really matter since there is 0 late coefficient.
                  return {{{}, {}}};
                }
                a64::vector<double> coeffs;
                coeffs.resize(nLateCoeffs);

                for(auto & c : async_convs) {
                    c->setup(
                    {
                        submission_period,
                        queue_size,
                        {
                            nAsyncScalesDropped
                        }
                    });
                    c->setCoefficients(coeffs);
                }
               if(phasing.mode == SimulationPhasingMode::On)
               {
                    int n=0;
                    int const total_sz = phasing.groupSize ? *phasing.groupSize : async_convs.size();
                    Assert(total_sz);
                    SetupParam sp{
                        1,
                        1,
                        {
                            CountDroppedScales(1)
                        }
                    };
                    for(auto & c : async_convs) {
                        dephase(total_sz,
                                n % total_sz,
                                sp,
                                *c);
                        ++n;
                    }
                }
      
                AudioHostSimulator simu{
                    frame_rate,
                    n_audio_frames_per_cb,
                    n_channels,
                    n_audio_channels
                };
                
                std::vector<QueueMetrics> m;
                m.reserve(async_convs.size());
                for(auto const & c : async_convs) {
                    int const metric_period = c->getAsyncAlgo().getBiggestScale();
                    m.emplace_back(metric_period,
                                   numPeriods);
                }
                
                auto p = simu.simulate([&async_convs, &m]
                         (std::vector<a64::vector<float>> const & inputs,
                          std::vector<a64::vector<float>> & outputs,
                          int const cur_frame,
                          int const n_frames) {
                    Assert(inputs.size() == async_convs.size());
                    Assert(!outputs.empty());
                    
                    for(int i=0; i < async_convs.size(); ++i) {
                        auto & conv = *async_convs[i];
                        auto const & input = inputs[i];
                        auto & output = outputs[i%outputs.size()];

                        bool const assign = i < outputs.size();

                        if(assign) {
                            conv.stepAssignVectorized(input.data()+cur_frame,
                                                      output.data()+cur_frame,
                                                      n_frames);
                        }
                        else {
                            conv.stepAddVectorized(input.data()+cur_frame,
                                                   output.data()+cur_frame,
                                                   n_frames);
                        }
              
                        
                        if(!m[i].recordQueueSize(conv.getResultQueueSize(),
                                                 conv.getSignalQueueSize(),
                                                 conv.countErrorsWorkerTooSlow(),
                                                 n_frames)) {
                            // the number of periods has elapsed or the worker is too slow
                            return false;
                        }
                    }
                    return true;
                });
                
                if(!p) {
                    // missed deadline
                    return {};
                }
                return std::make_pair(*p,m);
            }
        };
      
      template<typename T, typename FFTTag>
      struct PartitionAlgo< ZeroLatencyScaledAsyncConvolution<T, FFTTag> > {
          using Convolution = ZeroLatencyScaledAsyncConvolution<T, FFTTag>;
          using SetupParam = typename Convolution::SetupParam;
          using PS = PartitionningSpec<SetupParam>;
          using PSpecs = PartitionningSpecs<SetupParam>;

          using LateHandler = typename Convolution::LateHandler;

          static PSpecs run(int n_channels,
                            int n_audio_channels,
                            int n_audio_frames_per_cb,
                            int total_response_size,
                            double frame_rate,
                            double max_avg_time_per_sample,
                            std::ostream & os,
                            SimulationPhasing const & phasing)
        {
                // There is a balance to find between the risk of audio dropouts due to:
                //   - A : the earlyhandler having too many coefficients to handle
                //   - B : the latehandler having too many (small) scales to handle
                //
                // in the choice of nAsyncScalesDropped, here we chose to care about 'B':
                //
                // We set nAsyncScalesDropped such that, in the async part, we can use scale convolution
                // on the range of ffts where it's more optimal than brute force convolution.
                //
                // But using nDroppedOptimalFor_Split_Bruteforce_Fft instead of 0 for nAsyncScalesDropped
                // incurs more latency in the latehandler, hence more coefficients to handle in the early handler.
                //
                // So this design choice may be changed in the future.
              auto const nAsyncScalesDropped = ScaleConvolution_::nDroppedOptimalFor_Split_Bruteforce_Fft;

                auto pSpecsLate = PartitionAlgo<LateHandler>::run(n_channels,
                                                                  n_audio_channels,
                                                                  n_audio_frames_per_cb,
                                                                  total_response_size,
                                                                  frame_rate,
                                                                  max_avg_time_per_sample,
                                                                  os,
                                                                  phasing,
                                                                  nAsyncScalesDropped).getWithSpread();
                PS ps;
              if(pSpecsLate.cost) {
                ps.cost = {{{},{}}, *pSpecsLate.cost};
                ps.cost->setCost(pSpecsLate.getCost());
              }
              return {ps,{}};
          }
      };
      template<typename T, typename FFTTag>
      struct PartitionAlgo< ZeroLatencyScaledAsyncConvolutionOptimized<T, FFTTag> > {
          using Convolution = ZeroLatencyScaledAsyncConvolutionOptimized<T, FFTTag>;
          using SetupParam = typename Convolution::SetupParam;
          using PS = PartitionningSpec<SetupParam>;
          using PSpecs = PartitionningSpecs<SetupParam>;

          using EarlyHandler = typename Convolution::EarlyHandler;
          using LateHandler = typename Convolution::LateHandler;

          using LateSetupParam = typename LateHandler::SetupParam;
      
      static PSpecs run(int n_channels,
                        int n_audio_channels,
                        int n_audio_frames_per_cb,
                        int total_response_size,
                        double frame_rate,
                        double max_avg_time_per_sample,
                        std::ostream & os,
                        SimulationPhasing const & phasing)
      {
          os << "Late handler optimization:" << std::endl;

          // There is a balance to find between the risk of audio dropouts due to:
          //   - A : the earlyhandler having too many coefficients to handle
          //   - B : the latehandler having too many (small) scales to handle
          //
          // in the choice of nAsyncScalesDropped, here we chose to care about 'B':
          //
          // We set nAsyncScalesDropped such that, in the async part, we can use scale convolution
          // on the range of ffts where it's more optimal than brute force convolution.
          //
          // But using nDroppedOptimalFor_Split_Bruteforce_Fft instead of 0 for nAsyncScalesDropped
          // incurs more latency in the latehandler, hence more coefficients to handle in the early handler.
          //
          // So this design choice may be changed in the future.
        auto const nAsyncScalesDropped = ScaleConvolution_::nDroppedOptimalFor_Split_Bruteforce_Fft;

          auto indent = std::make_unique<IndentingOStreambuf>(os);
          auto pSpecsLate = PartitionAlgo<LateHandler>::run(n_channels,
                                                            n_audio_channels,
                                                            n_audio_frames_per_cb,
                                                            total_response_size,
                                                            frame_rate,
                                                            max_avg_time_per_sample,
                                                            os,
                                                            phasing,
                                                            nAsyncScalesDropped).getWithSpread();
          indent.reset();
          
          PS ps;
          if(pSpecsLate.cost) {
              int nEarlyCoeffs = deduceEarlyHandlerCoeffs(*(pSpecsLate.cost));
              os << "Deduced early handler coeffs:" << nEarlyCoeffs << std::endl;
              os << "Early handler optimization:" << std::endl;
              
              auto indent = std::make_unique<IndentingOStreambuf>(os);
              auto pSpecsEarly = PartitionAlgo<EarlyHandler>::run(n_channels,
                                                                  n_audio_channels,
                                                                  n_audio_frames_per_cb,
                                                                  nEarlyCoeffs,
                                                                  frame_rate,
                                                                  max_avg_time_per_sample,
                                                                  os).getWithSpread();
              indent.reset();
              
              if(pSpecsEarly.cost) {
                  ps.cost = {*pSpecsEarly.cost, *pSpecsLate.cost};
                  ps.cost->setCost(pSpecsEarly.getCost() + pSpecsLate.getCost());
              }
          }
        return {ps,{}};
      }
      private:
      static int deduceEarlyHandlerCoeffs(LateSetupParam p) {
        return p.getImpliedLatency();
      }
      };

namespace SameSizeScales {
    /*
     Assuming that scales have resolutions of 1,2,4,8...
     for every scale, the number of blocks = { N(k), k in 1 .. S }
     Assuming that between scale i and (i+1) there are nOverlap*2^(i-1) blocks in common
     N = sum (n in 0 .. (S-1)) 2^n N(n+1) - scaleFadeSz::inSmallerUnits * sum (n in 1.. (S-1)) 2^(n-1)
     and since we want the number of blocks to be equal
     (hence the namespace name 'SameSizeScales'), we will solve this:
     N = NF * sum (n in 0 .. (S-1)) 2^n - scaleFadeSz::inSmallerUnits * sum (n in 0.. (S-2)) 2^n
     N = NF * (2^S - 1) - scaleFadeSz::inSmallerUnits * (2^(S-1)-1)
     NF = (N + scaleFadeSz::inSmallerUnits * (2^(S-1)-1)) / (2^S - 1)
     */
    static inline int get_scale_sz(int response_sz, int n_scales) {
      assert(response_sz>=0);
      int numerator = response_sz + static_cast<int>(scaleFadeSz::inSmallerUnits * (static_cast<int>(pow2(n_scales-1)) - 1));
      int denominator = static_cast<int>(pow2(n_scales)) - 1;
      int res = ceil(numerator / static_cast<double>(denominator));
      if(subSamplingAllowsEvenNumberOfCoefficients || n_scales <= 1) {
        return res;
      }
      if(res%2) {
        return res + 1;
      }
      return res;
    }
    static inline int get_max_response_sz(int n_scales, int sz_one_scale) {
      return
      sz_one_scale * (pow2(n_scales) - 1) -
      scaleFadeSz::inSmallerUnits * (static_cast<int>(pow2(n_scales-1)) - 1);
    }
    
    int constexpr getDelays(int scale_sz, int partition_sz) {
      // split = nFadeIn + latB - latA
      // scale_sz = nFadeIn + (1 + 2*(delay + latA)) - latA
      // scale_sz - 1 - nFadeIn = 2*delay + latA
      // delay = 0.5 * (scale_sz - 1 - nFadeIn - latA)
      // delay = 0.5 * (scale_sz - 1 - nFadeIn - (2*partition_sz-1))
      // delay = (0.5 * (scale_sz - nFadeIn)) - partition_sz
      assert(0 == (scale_sz-scaleFadeSz::inSmallerUnits) % 2);
      return ((scale_sz-scaleFadeSz::inSmallerUnits) / 2) - partition_sz;
    }
  };

  

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

  
  constexpr int minLatencyLateHandlerWhenEarlyHandlerIsDefaultOptimizedFIRFilter = ScaleConvolution_::latencyForDroppedConvolutions(ScaleConvolution_::nDroppedOptimalFor_Split_Bruteforce_Fft);

  template<typename C>
  int constexpr lateHandlerLatency(int partition_sz) {
    using LateHandler = typename EarlyestDeepest<typename C::LateHandler>::type;
    auto res = LateHandler::getLatencyForPartitionSize(partition_sz);
    return res;
  }
  
  template<typename C>
  int constexpr getLateHandlerMinLg2PartitionSz() {
    using LateHandler = typename EarlyestDeepest<typename C::LateHandler>::type;

    int partition_sz = 1;
    for(;;partition_sz *= 2) {
      if(LateHandler::getLatencyForPartitionSize(partition_sz) >= minLatencyLateHandlerWhenEarlyHandlerIsDefaultOptimizedFIRFilter) {
        break;
      }
    }
    
    return power_of_two_exponent(partition_sz);
  }
  

        // Subsampling can be used to diminish the resolution of the impulse response tail,
        // it makes the computations use less CPU cycles:
        enum class ResponseTailSubsampling {
          // response is used at full resolution everywhere (most CPU intensive):
          FullRes, // 0
          ScaleCount_1 = FullRes,
          // the beginning of the response is at full resolution, then half resolution:
          UpToHalfRes, // 1
          ScaleCount_2 = UpToHalfRes,
          // the beginning of the response is at full resolution, then half resolution, then quarter resolution:
          UpToQuarterRes, // 2
          ScaleCount_3 = UpToQuarterRes,
          // the beginning of the response is at full resolution, then half resolution, then quarter resolution, then heighth resolution:
          UpToHeighthRes, // 3
          ScaleCount_4 = UpToHeighthRes,
          // If in doubt, use this mode: the least number of scales will be used
          // if we can afford the induced computations (using an auto optimizing algorithm).
          HighestAffordableResolution, // 4
        };

      template<typename Convolution>
      range<int> getScaleCountRanges(ResponseTailSubsampling rts) {
          range<int> r;
          
          if constexpr (Convolution::has_subsampling) {
              switch(rts) {
                  case ResponseTailSubsampling::ScaleCount_1:
                      r.extend(1);
                      break;
                  case ResponseTailSubsampling::ScaleCount_2:
                      r.extend(2);
                      break;
                  case ResponseTailSubsampling::ScaleCount_3:
                      r.extend(3);
                      break;
                  case ResponseTailSubsampling::ScaleCount_4:
                      r.extend(4);
                      break;
                  default:
                      assert(0);
                  case ResponseTailSubsampling::HighestAffordableResolution:
                      r.set(1, nMaxScales);
                      break;
              }
          }
          else {
              r.extend(1);
          }
          
          return r;
      }
      

      template<typename T, typename FFTTag, typename SP>
      auto mkSubsamplingSetupParams(SP const & params,
                                    int const n_scales,
                                    int const scale_sz) {
        using SetupParam = typename ZeroLatencyScaledFineGrainedPartitionnedSubsampledConvolution<T,FFTTag>::SetupParam;
        int delay = 0;
        if(n_scales > 1) {
          delay = SameSizeScales::getDelays(scale_sz, params.partition_size);
          Assert(delay > 0);
        }

        // set the delays
        
        // we disable the unused scales by setting the partition size to 0.
        auto zero = SP::makeInactive();
        static_assert(4==nMaxScales);
        std::array<SP, nMaxScales> ps {
          params,
          (n_scales >= 2)?params:zero,
          (n_scales >= 3)?params:zero,
          (n_scales >= 4)?params:zero
        };
        Assert(n_scales <= nMaxScales);
        return SetupParam
        {
            {{},{}},
          {
            ps[0],
            {
              (n_scales >= 2)?delay:0,
              {
                ps[1],
                {
                  (n_scales >= 3)?delay:0,
                  {
                    ps[2],
                    {
                      (n_scales >= 4)?delay:0,
                      ps[3]
                    }
                  }
                }
              }
            }
          }
        };
      }

      
  template<typename T, typename FFTTag>
  struct PartitionAlgo< ZeroLatencyScaledFineGrainedPartitionnedConvolution<T,FFTTag> > {
    using NonAtomicConvolution = ZeroLatencyScaledFineGrainedPartitionnedConvolution<T,FFTTag>;
    
  private:
    using LateHandler = typename NonAtomicConvolution::LateHandler;
  public:
      using SetupParam = typename NonAtomicConvolution::SetupParam;
    using PS = PartitionningSpec<SetupParam>;
    using PSpecs = PartitionningSpecs<SetupParam>;

    static PSpecs run(int n_channels,
                      int n_audio_channels,
                      int n_audio_frames_per_cb,
                      int total_response_size,
                      double frame_rate,
                      double max_avg_time_per_sample,
                      std::ostream & os) {
      if(total_response_size <= minLatencyLateHandlerWhenEarlyHandlerIsDefaultOptimizedFIRFilter) {
        // in that case there is no late handling at all.
      PS ps;
      ps.cost =
      {
        {{},{}},
        FinegrainedSetupParam::makeInactive()
      };
      ps.cost->setCost(0.f);
        return {ps,ps};
      }
      auto late_handler_response_size_for_partition_size =
      [total_response_size](int partition_size) -> Optional<int> {
        if(LateHandler::getLatencyForPartitionSize(partition_size) < minLatencyLateHandlerWhenEarlyHandlerIsDefaultOptimizedFIRFilter) {
          // invalid case. We pass 'getMinLg2PartitionSz()' to 'run' so that
          // this case doesn't happen on the first try.
          return {};
        }
        // we substract the count of coefficients handled by the early coefficients handler.
        // it favors long partitions because we don't take into account
        // the cost of the early coefficient handler, but for long responses, where optimization matters,
        // the induced bias is negligible.
        int n_coeffs_early_handler = lateHandlerLatency<NonAtomicConvolution>(partition_size);
        assert(n_coeffs_early_handler >= minLatencyLateHandlerWhenEarlyHandlerIsDefaultOptimizedFIRFilter);

        auto late_response_sz = std::max(0,total_response_size - n_coeffs_early_handler);

        if(late_response_sz <= 0) {
      // should we return {0} ?
          return {};
        }
        return {late_response_sz};
      };
      auto lateRes = PartitionAlgo<LateHandler>::run(n_channels,
                                             1,
                                             n_audio_frames_per_cb,
                                             late_handler_response_size_for_partition_size,
                                             getLateHandlerMinLg2PartitionSz<NonAtomicConvolution>(),
                                             os);
      PS withSpread;
      if(lateRes.with_spread.cost) {
          withSpread.cost = SetupParam
          {
              {{},{}},
              lateRes.with_spread.cost.value()
          };
          withSpread.cost->setCost(lateRes.with_spread.cost->getCost());
      }
      PS withoutSpread;
      if(lateRes.without_spread.cost) {
          withoutSpread.cost = SetupParam
          {
              {{},{}},
              lateRes.without_spread.cost.value()
          };
          withoutSpread.cost->setCost(lateRes.without_spread.cost->getCost());
      }
      return {
        withSpread,
        withoutSpread
      };
    }
  };
  
      template<typename T>
      struct WithNScales {
        T val;
        int n_scales;
      
        auto getCost() const { return val.getCost(); }
      
      void logReport(int n_channels,
                     double theoretical_max_avg_time_per_frame,
                     std::ostream & os) {
          val.logReport(n_channels, theoretical_max_avg_time_per_frame, os);
          os << "- using ";
          // TODO range of subsampling regions: highest quality region: xxx samples / 2-subsampled region : xxx samples / 4-subsampled
          if(n_scales <= 1) {
            os << "full tail resolution";
          }
          else {
            os << "reduced tail resolution with " << n_scales - 1 << " subsampling regions";
          }
      }
      };
      
      template<typename T, typename FFTTag>
      struct PartitionAlgo< ZeroLatencyScaledFineGrainedPartitionnedSubsampledConvolution<T,FFTTag> > {
        using NonAtomicConvolution = ZeroLatencyScaledFineGrainedPartitionnedSubsampledConvolution<T,FFTTag>;
        
      private:
        using LateHandler = typename EarlyestDeepest<typename NonAtomicConvolution::LateHandler>::type;
      public:
        using SetupParam = typename NonAtomicConvolution::SetupParam;//LateHandler::SetupParam;
        using PS = PartitionningSpec<SetupParam>;
        using PSpecs = PartitionningSpecs2<WithNScales<SetupParam>>;

      // il faudrait que a chaque fois ca renvoie PartitionningSpecs<NonAtomicConvolution::SetupParam> pour bien generaliser
        static PSpecs run(int const n_response_channels,
                          int const n_audio_channels,
                          int const n_audio_frames_per_cb,
                          int const n_response_frames,
                          double const frame_rate,
                          double const max_avg_time_per_sample,
                          std::ostream & os,
                          ResponseTailSubsampling rts) {
      
      Optional<PSpecs> res;
      
      range<int> const scales = getScaleCountRanges<NonAtomicConvolution>(rts);

      for(int n_scales = scales.getMin(); n_scales <= scales.getMax(); ++n_scales) {
        
        auto partit = runForScale(n_response_channels,
                                         n_audio_channels,
                                         n_audio_frames_per_cb,
                                         n_response_frames,
                                  n_scales,
                                         frame_rate,
                                         os);
        auto & part = partit.getWithSpread();
        if(!part.cost) {
          os << "Discard n_scales " << n_scales << std::endl;
          continue;
        }
        
        if(!res || (res->getCost() && (*res->getCost() > part.getCost()))) {
          res = {part, part};
        }

        if(!res || !res->getCost()) {
            continue;
        }
        if(*res->getCost() < max_avg_time_per_sample) {
          os << "Optimization criteria met with " << n_scales << " scaling levels." << std::endl;
          return *res;
        }
        os << "cost " << *res->getCost() << " >= " << max_avg_time_per_sample << std::endl;
      }
      throw std::runtime_error("Optimization criteria not met, there are not enough scaling levels.");

      }
  private:
      static PSpecs runForScale(int n_channels,
                   int n_audio_channels,
                   int n_audio_frames_per_cb,
                   int total_response_size,
                   int n_scales,
                   double frame_rate,
                   std::ostream & os)
      {
          if(total_response_size <= minLatencyLateHandlerWhenEarlyHandlerIsDefaultOptimizedFIRFilter) {
            // in that case there is no late handling at all.
            if(n_scales > 1) {
              LG(WARN, "not enough coefficients to use %d scales", n_scales);
              return {};
            }
            WithNScales<SetupParam> o
            {
                mkSubsamplingSetupParams<T, FFTTag>(FinegrainedSetupParam::makeInactive(), 1, 0),
                n_scales
            };
            o.val.setCost(0.f);
            return {{o},{o}};
          }
          auto scale_size_for_partition_size =
          [total_response_size, n_scales](int partition_size) -> Optional<int> {
            if(LateHandler::getLatencyForPartitionSize(partition_size) < minLatencyLateHandlerWhenEarlyHandlerIsDefaultOptimizedFIRFilter) {
              // invalid case. We pass 'getMinLg2PartitionSz()' to 'run' so that
              // this case doesn't happen on the first try.
              return {};
            }
            // we substract the count of coefficients handled by the early coefficients handler.
            // it favors long partitions because we don't take into account
            // the cost of the early coefficient handler, but for long responses, where optimization matters,
            // the induced bias is negligible.
            int n_coeffs_early_handler = lateHandlerLatency<NonAtomicConvolution>(partition_size);
            assert(n_coeffs_early_handler >= minLatencyLateHandlerWhenEarlyHandlerIsDefaultOptimizedFIRFilter);

            auto late_response_sz = std::max(0,total_response_size - n_coeffs_early_handler);

            auto scale_sz = SameSizeScales::get_scale_sz(late_response_sz, n_scales);
            if(scale_sz <= 0) {
              return {};
            }
            // TODO return also the size of early coefficients, and take the early cost into account.
            if(n_scales > 1) {
              // verify that we can have more than one scale (i.e the delay needed to scale is strictly positive)
              // in theory we could relax the constraint (0 delays are ok but the implementation
              // doesn't support that).
              if(SameSizeScales::getDelays(scale_sz, partition_size) <= 0) {
                return {};
              }
            }
            return {scale_sz};
          };
          auto lateRes = PartitionAlgo<LateHandler>::run(n_channels,
                                                 n_scales,
                                                 n_audio_frames_per_cb,
                                                 scale_size_for_partition_size,
                                                 getLateHandlerMinLg2PartitionSz<NonAtomicConvolution>(),
                                                 os);
          std::optional<WithNScales<SetupParam>> withSpread;
          if(lateRes.with_spread.cost) {
            std::optional<int> const mayScaleSz = scale_size_for_partition_size(lateRes.with_spread.cost->partition_size);
            if(mayScaleSz) {
                withSpread = {
                    mkSubsamplingSetupParams<T, FFTTag>(lateRes.with_spread.cost.value(), n_scales, mayScaleSz.value()),
                    n_scales
                };
                withSpread->val.setCost(lateRes.with_spread.cost->getCost());
            }
            else {
                throw std::logic_error("no scale size");
            }
          }
          std::optional<WithNScales<SetupParam>> withoutSpread;
          if(lateRes.without_spread.cost) {
            std::optional<int> const mayScaleSz = scale_size_for_partition_size(lateRes.without_spread.cost->partition_size);
            if(mayScaleSz) {
                withoutSpread = {
                  mkSubsamplingSetupParams<T, FFTTag>(lateRes.without_spread.cost.value(), n_scales, mayScaleSz.value()),
                  n_scales
                };
                withoutSpread->val.setCost(lateRes.without_spread.cost->getCost());
            }
          }
          return {
              {withSpread},
              {withoutSpread}
          };
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
