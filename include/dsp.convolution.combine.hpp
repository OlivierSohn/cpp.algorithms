
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

        void logSubReport(std::ostream & os) override {
            os << "a : " << std::endl;
            aParams.logSubReport(os);
            os << "b : " << std::endl;
            bParams.logSubReport(os);
        }
    };

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
   * The default value is optimal for the system I developped on (OSX / Macbook Air 2015 / intel core i7).
   * see 'TEST(BenchmarkConvolutions, scaled_vs_brute)' on how to determine the best value for a given system.
   *
   * TODO ideally this value should be a global, computed at initialization time.
   */
  constexpr int scaleConvolutionOptimalNDropped = 6;

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

    static constexpr int nCoefficientsFadeIn = 0;
    static_assert(A::nCoefficientsFadeIn == 0); // else we need to handle them

    // note that it wouldn't make much sense for 'A' to have subsampling in a ScaleConvolution
      static constexpr bool has_subsampling = A::has_subsampling;
      
    struct SetupParam : public Cost {
        SetupParam(int nDropped = scaleConvolutionOptimalNDropped)
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
      int nDropped;
      
      void logSubReport(std::ostream & os) override {
      os << "dropped: " << nDropped << std::endl;
      }
    };

    bool isZero() const {
      return v.empty();
    }
    void reset() {
      v.clear();
      x.clear();
      progress = 0;
      nDroppedConvolutions = 0;
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

    void applySetup(SetupParam const & p) {
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
          v.clear();
          return;
        }
      }
      x.resize(pow2(n)); // including padding for biggest convolution
    }

    bool isValid() const {
        if(nDroppedConvolutions < 0) {
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
        return pow2(nDroppedConvolutions)-1;
    }
      
      int countCoefficients() const {
          return pow2(nDroppedConvolutions + v.size()) - pow2(nDroppedConvolutions);
      }
      
      int getBiggestScale() const {
          if(v.empty()) {
              return 0;
          }
          return pow2(nDroppedConvolutions + v.size() - 1);
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
  
  template<typename T, typename FFTTag = fft::Fastest, LatencySemantic Lat = LatencySemantic::DiracPeak>
  using ZeroLatencyScaledFineGrainedPartitionnedConvolution =
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
      using ZeroLatencyScaledAsyncConvolution =
      SplitConvolution<
      // handle the first coefficients using a zero-latency filter:
        OptimizedFIRFilter<T, FFTTag>,
      // handle subsequent coefficients using an asynchronous scaling convolution
        AsyncCPUConvolution <
          ScaleConvolution <
            FFTConvolutionCore<T, FFTTag>
          >
        >
      >
      ;
      

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
                            int n_scales,
                            double frame_rate,
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

      template<typename T, typename FFTTag>
      struct PartitionAlgo< ZeroLatencyScaledAsyncConvolution<T, FFTTag> > {
          using Convolution = ZeroLatencyScaledAsyncConvolution<T, FFTTag>;
          using SetupParam = typename Convolution::SetupParam;
          using PS = PartitionningSpec<SetupParam>;
          using PSpecs = PartitionningSpecs<SetupParam>;

          struct Periodically {
              Periodically(std::chrono::steady_clock::duration time_step,
                           int n_max_periods = 1000000)
              : time_step(time_step)
              {
                  time_calls.reserve(n_max_periods);
                  t_zero = std::chrono::steady_clock::now();
              }
              
              template<typename F>
              void exec(F f) {
                  while(true) {
                      nextPeriod();
                      if(!f()) {
                          return;
                      }
                  }
              }
              
              std::chrono::steady_clock::time_point t_zero;
              std::chrono::steady_clock::duration time_step;
              std::vector<std::chrono::steady_clock::time_point> time_calls;
          private:
              void nextPeriod() {
                  auto t = time_step * i + t_zero;
                  ++i;
                  auto now = std::chrono::steady_clock::now();
                  if(now < t) {
                      auto dt = t - now;
                      std::this_thread::sleep_for(dt);
                  }
                  time_calls.push_back(std::chrono::steady_clock::now());
              };
              
              int i=0;
          };

          struct AudioHostSimulator {
              double const frame_rate;
              int const n_audio_frames_per_cb;
              int const n_audio_channels;
              
              // returns the size of the variation of the queue size
              template<typename F>
              Periodically simulate(F f) {
                  std::vector<std::vector<float>> inputs, outputs;
                  inputs.resize(n_audio_channels);
                  outputs.resize(n_audio_channels);
                  
                  for(auto & i : inputs) {
                      i.resize(n_audio_frames_per_cb);
                  }
                  for(auto & o : outputs) {
                      o.resize(n_audio_frames_per_cb);
                  }
                  
                  auto f2 = [&inputs, &outputs, f](){ return f(inputs, outputs); };

                  Periodically p(std::chrono::nanoseconds(static_cast<int>(1e9 * n_audio_frames_per_cb / frame_rate)));
                  p.exec(f2);
                  return p;
              }
          };
          
          struct Metrics {
              range<int> queueSizeRange;
              
              // if gradients are positive it means the worker is too slow:
              // the queue, on average, keeps increasing in size.

              double maxQueueSize_gradient = 0; // diff between first and last period over number of periods
              int maxQueueSize_max_local_gradient = 0; // diff between consecutive periods
              
              void mergeWorst(Metrics const & o) {
                  queueSizeRange.extend(o.queueSizeRange);
                  maxQueueSize_gradient = std::max(maxQueueSize_gradient,
                                                   o.maxQueueSize_gradient);
                  maxQueueSize_max_local_gradient = std::max(maxQueueSize_max_local_gradient,
                                                             o.maxQueueSize_max_local_gradient);
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
                  queueSizeRangeByPeriod.reserve(nPeriods);
                  queueSizeRangeByPeriod.resize(1);
              }
              
          private:
              int skip = 5;
              int const nPeriods;
              int const period_frames;
              int period = 0;
              int n_cur_frame = 0;
              
              std::vector<range<int>> queueSizeRangeByPeriod;
          public:
              bool recordQueueSize(int queueSize, int nFrames) {
                  assert(period < nPeriods);
                  if(skip) {
                      --skip;
                      return true;
                  }
                  n_cur_frame += nFrames;
                  if(n_cur_frame >= period_frames) {
                      n_cur_frame -= period_frames;
                      ++period;
                      if(period >= nPeriods) {
                          return false;
                      }
                      queueSizeRangeByPeriod.emplace_back();
                  }
                  queueSizeRangeByPeriod.back().extend(queueSize);
                  return true;
              }
              
              Metrics getMetrics() const {
                  Metrics m;
                  for(auto const & r : queueSizeRangeByPeriod) {
                      m.queueSizeRange.extend(r);
                  }
                  {
                      std::optional<int> prevMaxQueueSize;
                      std::optional<int> maxGradient;
                      for(auto const & r : queueSizeRangeByPeriod) {
                          if(!prevMaxQueueSize) {
                              prevMaxQueueSize = r.getMax();
                          }
                          else {
                              auto gradient = std::abs(*prevMaxQueueSize - r.getMax());
                              if(!maxGradient || *maxGradient < gradient) {
                                  maxGradient = gradient;
                              }
                          }
                      }
                      m.maxQueueSize_max_local_gradient = *maxGradient;
                  }
                  {
                      auto & first = queueSizeRangeByPeriod.front();
                      auto & last = queueSizeRangeByPeriod.back();
                      m.maxQueueSize_gradient = std::abs(first.getMax() - last.getMax()) / static_cast<double>(nPeriods);
                  }
                  return m;
              }
          };

          static PSpecs run(int n_channels,
                            int n_audio_channels,
                            int n_audio_frames_per_cb,
                            int total_response_size,
                            int n_scales,
                            double frame_rate,
                            std::ostream & os) {
              Assert(n_scales <= 1);
              if(n_scales > 1) {
                  throw std::logic_error("ZeroLatencyScaledAsyncConvolution doesn't support subsampling");
              }

              // optimize the variables of ZeroLatencyScaledAsyncConvolution, i.e:
              // - submission period
              // - queue size
              
              auto [periodically, queueMetrics] = computeQueueMetrics(n_channels,
                                                                      n_audio_channels,
                                                                      n_audio_frames_per_cb,
                                                                      total_response_size,
                                                                      frame_rate);
              std::vector<Metrics> metrics;
              metrics.reserve(queueMetrics.size());
              for(auto const & m : queueMetrics) {
                  metrics.push_back(m.getMetrics());
              }
              
              Metrics worst;
              for(auto const & m : metrics) {
                  worst.mergeWorst(m);
              }
              
              if(worst.maxQueueSize_gradient > 0.1) {
                  os << "Background worker is too slow";
                  return {};
              }
              
              int minNumFramesPerQueue = 4 * (1 + worst.queueSizeRange.getSpan());
              
              /*
               Pour masquer la latence due a l'asynchronicite il faut qu'il y aie un overlap de scales
               entre les scalings asynchrones et synchrones:

               synchrone
                     asynchrone
               
                     ...
                     8192
                     4096
                     2048
               1024  1024
               512   512
               256
               128
               64
               32
               16
               8
               4
               2
               1
               
               Si le nombre de frames par callback est faible, il se peut que l'on doive utiliser fingrainedpartition
               pour lisser les pics (et du coup changer le type de reverb).
               */
              

              // ideally, all samples should be submitted at the end of the callback.
              int submission_period = std::max(1,
                                               static_cast<int>(floor_power_of_two(n_audio_frames_per_cb)/4));
              int latency = 10000;
              int queue_size = latency/submission_period;
              
              int const nAsyncScaledDropped = 9;
              int const minAsynchronousFftSize = pow2(nAsyncScaledDropped);
              
              PS minimalPs;
              minimalPs.cost = SetupParam({
                  {},{}
                  
              },{
                  submission_period,
                  queue_size,
                  {
                      nAsyncScaledDropped
                  }
              });
              minimalPs.cost->setCost(0.);
              return {
                  minimalPs,
                  minimalPs
              };
          }
      private:
          static std::pair<Periodically,std::vector<QueueMetrics>>
          computeQueueMetrics(int n_channels,
                              int n_audio_channels,
                              int n_audio_frames_per_cb,
                              int total_response_size,
                              double frame_rate)
          {
              using AsyncPart = typename Convolution::LateHandler;
              std::vector<std::unique_ptr<AsyncPart>> async_convs;
              async_convs.reserve(n_channels);
              for(int i=0; i<n_channels; ++i) {
                  async_convs.push_back(std::make_unique<AsyncPart>());
              }

              // todo dephase

              int const num_frames_in_queue = 2*total_response_size; // total_response_size should be enough

              int const submission_period = n_audio_frames_per_cb;
              int queue_size = num_frames_in_queue / submission_period;
              if(queue_size * submission_period < num_frames_in_queue) {
                  ++queue_size;
              }
              
              // this is pessimistic, we can adjust this later based on the results of the simulation
              int const min_num_scales_dropped = power_of_two_exponent(n_audio_frames_per_cb);

              int const nDroppedCoeffs = std::max(0,
                                                  static_cast<int>(pow2(min_num_scales_dropped)) - 1);
              int const nLateCoeffs = total_response_size-nDroppedCoeffs;
              a64::vector<double> coeffs;
              coeffs.resize(nLateCoeffs);

              for(auto & c : async_convs) {
                  applySetup(*c,
                  {
                      submission_period,
                      queue_size,
                      {
                          min_num_scales_dropped
                      }
                  });
                  c->setCoefficients(coeffs);
              }
              AudioHostSimulator simu{
                  frame_rate,
                  n_audio_frames_per_cb,
                  n_audio_channels
              };
              
              std::vector<QueueMetrics> m;
              m.reserve(async_convs.size());
              for(auto const & c : async_convs) {
                  int const metric_period = c->getAsyncAlgo().getBiggestScale();
                  m.emplace_back(metric_period,
                                 4);
              }
              
              auto p = simu.simulate([&async_convs, &m]
                       (std::vector<std::vector<float>> const & inputs,
                        std::vector<std::vector<float>> & outputs) {
                  assert(inputs.size() == outputs.size());
                  assert(inputs.size());
                  int const n_frames = inputs[0].size();
                  int const n = async_convs.size() / inputs.size();
                  // n = 2 for true stereo
                  
                  for(int i=0; i < async_convs.size(); ++i) {
                      auto & conv = *async_convs[i];
                      auto const & input = inputs[i/n];
                      auto & output = outputs[i%n];

                      bool const assign = i < inputs.size();

                      if(assign) {
                          conv.stepAssignVectorized(input.data(),
                                                    output.data(),
                                                    input.size());
                      }
                      else {
                          conv.stepAddVectorized(input.data(),
                                                 output.data(),
                                                 input.size());
                      }

                      if(!m[i].recordQueueSize(conv.getResultQueueSize(),
                                               n_frames)) {
                          // this means the end of the simulation
                          return false;
                      }
                  }
                  return true;
              });
              
              return std::make_pair(p,m);
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

  
  constexpr int minLatencyLateHandlerWhenEarlyHandlerIsDefaultOptimizedFIRFilter = static_cast<int>(pow2(scaleConvolutionOptimalNDropped))-1;

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
  
  template<typename T, typename FFTTag>
  struct PartitionAlgo< ZeroLatencyScaledFineGrainedPartitionnedConvolution<T,FFTTag> > {
    using NonAtomicConvolution = ZeroLatencyScaledFineGrainedPartitionnedConvolution<T,FFTTag>;
    
  private:
    using LateHandler = typename EarlyestDeepest<typename NonAtomicConvolution::LateHandler>::type;
  public:
    using SetupParam = typename LateHandler::SetupParam;
    using PS = PartitionningSpec<SetupParam>;
    using PSpecs = PartitionningSpecs<SetupParam>;

    static PSpecs run(int n_channels,
                      int n_audio_channels,
                      int n_audio_frames_per_cb,
                      int total_response_size,
                      int n_scales,
                      double frame_rate,
                      std::ostream & os) {
      if(total_response_size <= minLatencyLateHandlerWhenEarlyHandlerIsDefaultOptimizedFIRFilter) {
        // in that case there is no late handling at all.
        if(n_scales > 1) {
          LG(WARN, "not enough coefficients to use %d scales", n_scales);
          return {};
        }
        PS ps;
        ps.cost = FinegrainedSetupParam::makeInactive();
        return {ps,ps};
      }
      auto late_handler_response_size_for_partition_size =
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
      return PartitionAlgo<LateHandler>::run(n_channels,
                                             n_scales,
                                             n_audio_frames_per_cb,
                                             late_handler_response_size_for_partition_size,
                                             getLateHandlerMinLg2PartitionSz<NonAtomicConvolution>(),
                                             os);
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
