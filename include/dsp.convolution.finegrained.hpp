
namespace imajuscule
{
  /*
   * In dsp.convolution.hpp we see this algorithm:
   *
   * The impulse response h is split in parts of equal length h1, h2, ... hn
   *
   *                     FFT(h1)
   *                        v
   *       +-----+    +-----------+  +-----+  +------+
   * x +-->| FFT |-|->| cplx mult |->| Add |->| IFFT |--> y
   *       +-----+ v  +-----------+  +-----+  +------+
   *          +--------+               ^
   *          |Delay(N)| FFT(h2)       |
   *          +--------+    v          |
   *               |  +-----------+    |
   *               |->| cplx mult |----|
   *               |  +-----------+    |
   *               .         .         .
   *               .         .         .
   *               .                   .
   *               v                   |
   *          +--------+               |
   *          |Delay(N)| FFT(hn)       |
   *          +----+---+    v          |
   *               |  +-----------+    |
   *               -->| cplx mult |----|
   *                  +-----------+
   *
   * 'PartitionnedFFTConvolutionCRTP' carries this monolithic computation
   * when needed, with no buffering, with a latency of the size of a partition.
   *
   * When partition sizes are an order of magnitude bigger than the audio callback buffer
   * (it is a very common case, since large partition
   * sizes is what makes this algorithm efficient on average),
   * many successive callback calls have very little work to do, and suddenly
   * a single callback call has to perform this huge monolithic computation.
   *
   * This is problematic because this audio callback can miss its deadline,
   * and then we'll hear a loud audio crack.
   *
   * To fix this, at the cost of a longer latency (twice the size of a partition size),
   * here we split the computation in several "grains" that can be computed at different times:
   *
   * - First there is the "FFT" grain which computes the fft of a chunk of the input signal
   *     (and does a little more than that, see the code)
   * - Then there are "multiplication" grains, where we multiply some delayed FFT of the input signal
   *    by some of the the FFT of the impulse response.
   *    The (max) number of vector multiplications per grain is the "multiplication group size".
   * - Finally there is the "IFFT" grain where we sum the results of the multiplications
   *    and do its inverse fft.
   *
   * An optimization algorithm minimizes the worst case cost for a single callback,
   * based on :
   *  - The callback buffer size (which we assume will not change in the future)
   *  - The count of simultaneous convolutions happening in a callback
   * the optimization algorithm gives the optimal parameters:
   *  - The size of the partitions
   *  - The count of multiplications per multiplication grain
   *  - The "phasing" of simultaneous convolutions to best interleave
   *      high-cost grains.
  */

  enum class GrainType {
    FFT,
    IFFT,
    MultiplicationGroup,
    Nothing
  };

  struct GrainsCosts {
    float fft, ifft, mult;
      
      bool operator == (GrainsCosts const & o) const {
          return fft == o.fft &&
          ifft == o.ifft &&
          mult == o.mult;
      }
  };


  struct FinegrainedSetupParam : public Cost {
    FinegrainedSetupParam() {}

    FinegrainedSetupParam(int partitionSz, int multiplication_group_size, int phase) : Cost(phase),
    multiplication_group_size(multiplication_group_size),
    partition_size(partitionSz)
    {}

    void setGrainsCosts(GrainsCosts gcosts) { grains_costs = gcosts; }

    int multiplication_group_size = 0;
    int partition_size = 0;
    GrainsCosts grains_costs;
      
      void logSubReport(std::ostream & os) const override {
          os << "Finegrained, fft size: " << partition_size << std::endl;
      }


    static FinegrainedSetupParam makeInactive() {
      FinegrainedSetupParam res{0,0,0};
      res.setCost(0.f);
      return res;
    }
      
      bool operator ==(FinegrainedSetupParam const & o) const {
          return o.multiplication_group_size == multiplication_group_size &&
          o.partition_size == partition_size &&
          o.grains_costs == grains_costs &&
          getCost() == o.getCost();
      }
  };



  template <typename Parent>
  struct FinegrainedFFTConvolutionBase : public Parent {
    friend class Test;

    using T = typename Parent::FPT;
    using FPT = T;
    using Tag = typename Parent::FFTTag;

    static constexpr int nComputePhaseable = 1;
    static constexpr int nCoefficientsFadeIn = 0;
      static constexpr bool has_subsampling = false;

    using SetupParam = FinegrainedSetupParam;

    using Parent::set_partition_size;
    using Parent::getLatencyForPartitionSize;
    using Parent::doLogComputeState;

      void logComputeState(std::ostream & os) const {
          os << "Finegrained ";
          if(isZero()) {
              os << "zero" << std::endl;
          }
          else {
              os << " grain " << grain_counter << "/" << countGrains()
              << " progress " << x.size() << "/" << getBlockSize() << std::endl;
              doLogComputeState(os);
          }
      }

    void setup(SetupParam const & p) {
      set_partition_size(p.partition_size);
      setMultiplicationGroupLength(p.multiplication_group_size);
    }
      std::array<int, 1> getComputePeriodicities() const {
          return {getBlockSize()};
      }
      // in [0, getComputePeriodicity())
      std::array<int, 1> getComputeProgresses() const {
          auto const sz = getBlockSize();
          auto const res = sz ? static_cast<int>(x.size() % sz) : 0;
          return {res};
      }
      void setComputeProgresses(std::array<int, 1> const & progresses) {
          if(!getBlockSize()) {
              return;
          }
          auto const p = progresses[0];
          while(getComputeProgresses()[0] != p) {
              step(0);
          }
      }

    static constexpr auto zero_signal = fft::RealSignal_<Tag, FPT>::zero;
    static constexpr auto copy = fft::RealSignal_<Tag, FPT>::copy;
    static constexpr auto get_signal = fft::RealSignal_<Tag, FPT>::get_signal;
    static constexpr auto add_scalar_multiply = fft::RealSignal_<Tag, FPT>::add_scalar_multiply;
    using RealSignal = typename fft::RealSignal_<Tag, FPT>::type;
    using Algo = typename fft::Algo_<Tag, FPT>;
    using Contexts = fft::Contexts_<Tag, FPT>;

    using Parent::get_fft_length;
    using Parent::getBlockSize;
    using Parent::doSetCoefficients;

    using Parent::getGranularMinPeriod;
    using Parent::doSetMultiplicationGroupLength;
    using Parent::isValid;
    using Parent::clear;
    using Parent::countPartitions;
    using Parent::countGrains;
    using Parent::getGrainNumber;
    using Parent::increment_grain;
    using Parent::compute_x_fft;
    using Parent::do_some_multiply_add;
    using Parent::get_multiply_add_result;

    FinegrainedFFTConvolutionBase() {
      it = y.end();
    }

    void setCoefficients(a64::vector<T> coeffs_) {
      if(coeffs_.size() < 2) {
        coeffs_.resize(2); // avoid ill-formed cases
      }
      auto const N = coeffs_.size();
      auto const fft_length = get_fft_length(N);
      fft.setContext(Contexts::getInstance().getBySize(fft_length));

      assert(fft_length > 0);

      result.resize(fft_length);
      x.reserve(fft_length);
      {
        y.resize(fft_length);
      }

      doSetCoefficients(fft, std::move(coeffs_));
      reset_states();
    }

    void reset() {
      y.clear();
      result.clear();
      reset_states();
      clear();
    }
    void flushToSilence() {
        Parent::reset_base_states();
        if(!result.empty()) {
            zero_signal(result);
        }
        x.clear();
        if(!y.empty()) {
            zero_signal(y);
        }
        it = y.begin();
        grain_counter = 0;
    }

    bool isZero() const {
      return result.empty();
    }

    void setMultiplicationGroupLength(int l) {
      doSetMultiplicationGroupLength(l);
      reset_states();
    }

  private:
    void reset_states() {
      Parent::reset_base_states();
      x.clear();
      if(!y.empty()) {
        zero_signal(y);
      }
      it = y.begin();
      grain_counter = 0;
    }

  public:

    // Just used for calibration
    void fastForwardToComputation(GrainType t, T val = 1) {
      switch(t) {
        case GrainType::FFT:
          while(x.size() != getBlockSize()-1) {
            step(val);
          }
          break;
        case GrainType::IFFT:
          while(getGrainNumber() != countGrains() - 1) {
            step(val);
          }
          while(grain_counter != getGranularMinPeriod() - 1) {
            step(val);
          }
          break;
        case GrainType::MultiplicationGroup:
          while(!x.empty()) {
            step(val);
          }
          while(grain_counter != getGranularMinPeriod() - 1) {
            step(val);
          }
          break;
      }
      assert(willComputeNextStep());
    }
    // Just used for calibration
    bool willComputeNextStep() const {
      return (x.size() == getBlockSize()-1) || (grain_counter+1 == getBlockSize()/countGrains());
    }
      
      T step(T val) {
          assert(isValid());
          auto g = nextGrain();
          assert(g.first > 0);

          x.emplace_back(val);
          ++it;
          ++grain_counter;
          
          if(g.first == 1 && g.second) {
            doGrain(*g.second);
          }
          
          assert(it < y.begin() + getBlockSize());
          assert(it >= y.begin());
          assert(!y.empty());
          
          return get_signal(*it);
      }
      
      template<typename FPT2, typename F>
      void vectorized(FPT2 const * input_buffer,
                      FPT2 * output_buffer,
                      int nSamples,
                      F f)
      {
          while(true) {
              assert(nSamples >= 0);
              auto g = nextGrain();
              if(nSamples < g.first) {
                  // we don't have enough samples to reach the grain
                  g.first = nSamples;
                  g.second.reset();
              }
              nSamples -= g.first;
              grain_counter += g.first;
              for(int i=0; i<g.first; ++i) {
                  x.emplace_back(input_buffer[i]);
              }
              for(int i=0; i<g.first; ++i) {
                  ++it;
                  if(i == g.first-1 && g.second) {
                      doGrain(*g.second);
                  }
                  f(output_buffer[i], get_signal(*it));
              }
              if(0 == nSamples) {
                  return;
              }
              input_buffer += g.first;
              output_buffer += g.first;
          }
      }
      
      template<typename FPT2>
      void stepAddVectorized(FPT2 const * input_buffer,
                             FPT2 * output_buffer,
                             int nSamples)
      {
          vectorized(input_buffer, output_buffer, nSamples,
                     [](FPT2 & output, FPT result){
              output += result;
          });
      }
      template<typename FPT2>
      void stepAssignVectorized(FPT2 const * input_buffer,
                             FPT2 * output_buffer,
                             int nSamples)
      {
          vectorized(input_buffer, output_buffer, nSamples,
                     [](FPT2 & output, FPT result){
              output = result;
          });
      }
      template<typename FPT2>
      void stepAddInputZeroVectorized(FPT2 * output_buffer,
                                      int nSamples)
      {
          while(true) {
              assert(nSamples >= 0);
              auto g = nextGrain();
              if(nSamples < g.first) {
                  // we don't have enough samples to reach the grain
                  g.first = nSamples;
                  g.second.reset();
              }
              nSamples -= g.first;
              grain_counter += g.first;
              for(int i=0; i<g.first; ++i) {
                  x.push_back({});
              }
              for(int i=0; i<g.first; ++i) {
                  ++it;
                  if(i == g.first-1 && g.second) {
                      doGrain(*g.second);
                  }
                  output_buffer[i] += get_signal(*it);
              }
              if(0 == nSamples) {
                  return;
              }
              output_buffer += g.first;
          }
      }

  private:
      
      std::pair<int, std::optional<GrainType>> nextGrain() const {
          int xSize = static_cast<int>(x.size());
          int const block_size = getBlockSize();
          int distanceToFFTGrain = block_size - xSize;
          assert(distanceToFFTGrain >= 0);
          
          int const n_grains = countGrains();
          int const granularity = block_size/n_grains;
          int distanceToOtherGrain = granularity - grain_counter;
          assert(distanceToOtherGrain >= 0);
          
          // in case of equality, FFT wins.
          if(distanceToOtherGrain < distanceToFFTGrain) {
              auto cur_grain = getGrainNumber();
              assert(cur_grain <= n_grains);
              if( cur_grain < n_grains - 1 ) {
                  return {distanceToOtherGrain, GrainType::MultiplicationGroup};
              }
              else if(cur_grain == n_grains - 1) {
                  return {distanceToOtherGrain, GrainType::IFFT};
              }
              else {
                  // spread not optimal
                  return {distanceToOtherGrain, GrainType::Nothing};
              }
          }
          else {
              return {distanceToFFTGrain, GrainType::FFT};
          }
      }
      
      void doGrain(GrainType g)
      {
          switch(g)
          {
              case GrainType::FFT:
              {
                  auto const block_size = getBlockSize();
                  assert(it == y.begin() + block_size); // make sure 'rythm is good', i.e we exhausted the first half of the y vector
                  it = y.begin();
                  
                  //////////////////////// FFT grain /////////////////////////////////////
                  
                  // pad x
                  
                  x.resize(get_fft_length());
                  
                  compute_x_fft(fft, x);
                  
                  x.clear();
                  
                  // in the same "grain of computation" we do the following.
                  // if it is too costly we must delay the computation
                  // in another grain and have a second y vector
                  
                  auto factor = 1 / (Algo::scale * Algo::scale * static_cast<T>(get_fft_length()));
                  
                  auto it_res = result.begin();
                  auto it_y = y.begin();
                  auto it_y_prev = it_y + block_size;
                  
                  // y = mix first part of result with second part of previous result
                  //
                  // 'first part of y' = factor * ('second part of y' + 'first part of result')
                  add_scalar_multiply(it_y, /* = */
                                      /* ( */ it_res, /* + */ it_y_prev /* ) */, /* x */ factor,
                                      block_size);
                  
                  // store second part of result for later
                  //
                  // 'second part of y' = 'second part of result'
                  copy(it_y   + block_size,
                       it_res + block_size,
                       block_size);
                  
                  increment_grain();
                  break;
              }
              case GrainType::MultiplicationGroup:
              {
                  do_some_multiply_add();
                  increment_grain();
                  break;
              }
              case GrainType::IFFT:
              {
                  fft.inverse(get_multiply_add_result(),
                              result,
                              get_fft_length());
                  increment_grain();
                  break;
              }
              case GrainType::Nothing:
                  break;
          }
          grain_counter = 0;
      }
      
  private:
    int grain_counter = 0;
    Algo fft;
    RealSignal x, y, result;
    typename decltype(y)::iterator it;
  };

  constexpr int countPartitions(int nCoeffs, int partition_size) {
    auto n_partitions = nCoeffs/partition_size;
    if(n_partitions * partition_size != nCoeffs) {
      assert(n_partitions * partition_size < nCoeffs);
      ++n_partitions;
    }
    return n_partitions;
  }


  template <typename T, typename Tag>
  struct FinegrainedPartitionnedFFTConvolutionCRTP {
    using FPT = T;
    using FFTTag = Tag;

    using RealSignal = typename fft::RealSignal_<Tag, FPT>::type;
    using Signal_value_type = typename fft::RealSignal_<Tag, FPT>::value_type;

    using CplxFreqs = typename fft::RealFBins_<Tag, FPT>::type;
    static constexpr auto zero = fft::RealFBins_<Tag, FPT>::zero;
    static constexpr auto multiply_add = fft::RealFBins_<Tag, FPT>::multiply_add;

    using Algo = typename fft::Algo_<Tag, FPT>;
      void doLogComputeState(std::ostream & os) const {
          os << countPartitions() << " partitions "
          << getMultiplicationsGroupMaxSize() << " multGroupMaxSize" << std::endl;
      }
      
    auto get_fft_length() const { return 2 * partition_size; }
    auto get_fft_length(int) const { return get_fft_length(); }

    bool empty() const {
      return ffts_of_partitionned_h.empty();
    }
    void clear() {
      ffts_of_partitionned_h.clear();
    }

    auto getBlockSize() const { return partition_size; }
    static constexpr int getLatencyForPartitionSize(int partSz) {
      return 2*partSz - 1;
    }
    auto getLatency() const { return getLatencyForPartitionSize(partition_size); }
    auto getGranularMinPeriod() const { return getBlockSize() / countGrains(); }
    bool isValid() const { return partition_size > 0 && mult_grp_len > 0 && countGrains() <= getBlockSize(); }

    int countPartitions() const { return ffts_of_partitionned_h.size(); }

    double getEpsilon() const {
      return countPartitions() * (fft::getFFTEpsilon<FPT>(get_fft_length()) + 2 * std::numeric_limits<FPT>::epsilon());
    }

  protected:
    void doSetMultiplicationGroupLength(int length) {
      mult_grp_len = length;
    }

    void reset_base_states() {
      grain_number = 1;
    }

    // 0 means that it will not be used.
    void set_partition_size(int sz) {
      partition_size = sz;
      assert(sz == 0 || is_power_of_two(sz));
    }

  public:
    auto getMultiplicationsGroupMaxSize() const { return mult_grp_len; }
    auto countMultiplicativeGrains() const { return 1 + (countPartitions()-1)/getMultiplicationsGroupMaxSize(); }
    static constexpr auto countNonMultiplicativeGrains() { return 2; }
    auto countGrains() const { return countNonMultiplicativeGrains() + countMultiplicativeGrains(); }

    int getLowestValidMultiplicationsGroupSize() const {
      // lowest valid mult_grp_len verifies:

      // countGrains() == getBlockSize()
      // 2 + 1 + (ffts_of_partitionned_h.size() - 1)/mult_grp_len == partition_size

      auto constexpr min_number_grains = countNonMultiplicativeGrains() + 1;
      auto diff = partition_size - min_number_grains;
      if(diff < 0) {
        // invalid configuration
        return getHighestValidMultiplicationsGroupSize();
      }
      assert(diff >=0 );
      for(int i=1;; ++i) {
        if( (static_cast<int>(ffts_of_partitionned_h.size()) - 1)/i <= diff) {
          return i;
        }
      }
    }

    int getHighestValidMultiplicationsGroupSize() const { return countPartitions(); }

    void doSetCoefficients(Algo const & fft, a64::vector<T> coeffs_) {
      assert(partition_size > 0);

      auto const n_partitions = imajuscule::countPartitions(coeffs_.size(), partition_size);
      // if one partition is partial, it will be padded with zeroes.
      coeffs_.resize(n_partitions * partition_size);

      ffts_of_delayed_x.resize(n_partitions);
      ffts_of_partitionned_h.resize(n_partitions);

      auto const fft_length = get_fft_length();

      for(auto & fft_of_partitionned_h : ffts_of_partitionned_h) {
        fft_of_partitionned_h.resize(fft_length);
      }
      for(auto & fft_of_delayed_x : ffts_of_delayed_x) {
        fft_of_delayed_x.resize(fft_length);
      }

      work.resize(fft_length);

      // compute fft of padded impulse response

      auto it_coeffs = coeffs_.begin();
      {
        RealSignal coeffs_slice(fft_length, Signal_value_type(0)); // initialize with zeros (second half is padding)
        for(auto & fft_of_partitionned_h : ffts_of_partitionned_h) {
          auto end_coeffs = it_coeffs + partition_size;
          assert(end_coeffs <= coeffs_.end());
          auto slice_it = coeffs_slice.begin();
          for(;it_coeffs != end_coeffs; ++it_coeffs, ++slice_it) {
            using RealT = typename RealSignal::value_type;
            *slice_it = RealT(*it_coeffs);
          }

          // coeffs_slice is padded with 0, because it is bigger than partition_size
          // and initialized with zeros.
          fft.forward(coeffs_slice.begin(), fft_of_partitionned_h, fft_length);
        }
      }
      assert(it_coeffs == coeffs_.end());
    }

  protected:

    int getGrainNumber() const { return grain_number; }

    void compute_x_fft(Algo const & fft, RealSignal const & x) {
      assert(grain_number == countGrains());
      grain_number = 0;
      auto & oldest_fft_of_delayed_x = *ffts_of_delayed_x.cycleEnd();
      auto const fft_length = get_fft_length();
      assert(fft_length == oldest_fft_of_delayed_x.size());
      fft.forward(x.begin(), oldest_fft_of_delayed_x, fft_length);
      ffts_of_delayed_x.advance();
    }

    void do_some_multiply_add() {
      auto const M = getMultiplicationsGroupMaxSize();
      auto const offset_base = M * (grain_number - 1);
      assert(offset_base >= 0);
      if(0 == offset_base) {
        zero(work);
      }

      assert(offset_base < ffts_of_partitionned_h.size());
      auto it_fft_of_partitionned_h = ffts_of_partitionned_h.begin() + offset_base;
      for(auto offset = offset_base, offset_end = std::min((offset_base+M), static_cast<int>(ffts_of_partitionned_h.size()));
          offset != offset_end;
          ++offset, ++it_fft_of_partitionned_h) {
        assert(it_fft_of_partitionned_h < ffts_of_partitionned_h.end());

        auto const & fft_of_delayed_x = ffts_of_delayed_x.get_backward(offset);

        multiply_add(work                /*   +=   */,
                     fft_of_delayed_x,   /*   x   */   *it_fft_of_partitionned_h);
      }
    }

    auto const & get_multiply_add_result() const {
      return work;
    }

    void increment_grain() {
      ++grain_number;
    }

  private:
    int mult_grp_len = 0;
    int partition_size = -1;
    int grain_number = 0;
    cyclic<CplxFreqs> ffts_of_delayed_x;
    std::vector<CplxFreqs> ffts_of_partitionned_h;

    CplxFreqs work;
  };

  template <typename T, typename FFTTag = fft::Fastest>
  using FinegrainedPartitionnedFFTConvolution = FinegrainedFFTConvolutionBase< FinegrainedPartitionnedFFTConvolutionCRTP<T, FFTTag> >;

  /*
   input parameters :
   - n_frames_audio_cb, n_channels, with_spread (those 3 can be reduced to 'equivalent_n_frames_cb')
   - impulse response length

   output parameters:
   - lg(partition size)

   1D - Gradient descent according to cost 'max grain computation time' with variable parameters 'lg(partition_size)'
   deducing 'number of multiplications per grain' by finding the parameter that leads to grain computations time just below max(fft, ifft),
   with the constraint that one computation at most occurs per 'equivalent' audio callback.
   */
  template<typename NonAtomicConvolution, typename SetupParam = typename NonAtomicConvolution::SetupParam, typename GradientDescent = GradientDescent<SetupParam>>
  void find_optimal_partition_size_for_nonatomic_convolution(GradientDescent & gradient_descent,
                                                               int n_iterations,
                                                               int n_channels,
                                                               int n_scales,
                                                               int n_frames,
                                                               std::function<Optional<int>(int)> n_coeffs_for_partition_sz,
                                                             int min_lg2_partitionsz, // must yield a valid result
                                                               bool constraint,
                                                               Optional<SetupParam> & min_val,
                                                             std::ostream & os) {
    //std::cout << "main thread: " << std::endl;
    //thread::logSchedParams();

    gradient_descent.setFunction( [n_frames, n_coeffs_for_partition_sz, n_scales, constraint, n_iterations, n_channels, &os] (int const lg2_partition_size, auto & val){
      using namespace profiling;
      using namespace std;
      using namespace std::chrono;

        // repeat 3 times (and take the min duration):
        // - first time used to cache things,
        // - then we do 2 others so that if one is preempted, the other is (likely) not
      constexpr auto n_atoms_repeat = 3;
      constexpr auto n_atoms_repeat_warmup = 0; // no need to warmup, because we take the _minimum_ duration

      if(lg2_partition_size < 0) {
        return ParamState::OutOfRange;
      }
      if(lg2_partition_size > 20) {
        throw logic_error("Gradient descent is not working?");
      }
      int const partition_size = pow2(lg2_partition_size);
      //            cout << "partition size : " << partition_size << endl;

      auto maybe_impulse_sz = n_coeffs_for_partition_sz(partition_size);
      if(!maybe_impulse_sz) {
        return ParamState::OutOfRange;
      }
      int const length_impulse = *maybe_impulse_sz;

      struct Test {
        using T = typename NonAtomicConvolution::FPT;
        Test(int partition_size, int length_impulse) {
          a64::vector<T> coeffs;
          coeffs.reserve(length_impulse);
          std::array<T,4> values {0.9,0.5,-0.2,-0.5};
          for(int i=0; i<length_impulse; ++i) {
            coeffs.push_back(values[i%values.size()]);
          }

          // the value for multiplication group size is not very important (it will be overriden later on)
          // but needs to lead to a valid convolution. We use the highest valid number:
          pfftcv.setup({partition_size, countPartitions(coeffs.size(), partition_size), 0});
          pfftcv.setCoefficients(std::move(coeffs));
        }

        bool isValid() const { return pfftcv.isValid(); }

        int getGranularMinPeriod() const { return pfftcv.getGranularMinPeriod(); }
        int countMultiplicativeGrains() const { return pfftcv.countMultiplicativeGrains(); }
        int getHighestValidMultiplicationsGroupSize() const { return pfftcv.getHighestValidMultiplicationsGroupSize(); }
        int getLowestValidMultiplicationsGroupSize() const { return pfftcv.getLowestValidMultiplicationsGroupSize(); }

        void setMultiplicationGroupLength(int mult_grp_length) {
          pfftcv.setMultiplicationGroupLength(mult_grp_length);
        }

        void prepare(GrainType g) {
          pfftcv.fastForwardToComputation(g);
        }
        void run() {
          assert(pfftcv.willComputeNextStep());
          pfftcv.step(1.f);
        }
      private:
        NonAtomicConvolution pfftcv;
      };

      // prepare tests

      // TODO [early coefficients cost] substract the early coefficients from length_impulse
      Test test(partition_size, length_impulse);
        if(!test.isValid()) {
          return ParamState::OutOfRange;
        }

      constexpr auto n_non_multiplicative_grains = NonAtomicConvolution::countNonMultiplicativeGrains();
      static_assert(2 == n_non_multiplicative_grains);
      static constexpr auto index_fft = 0;
      static constexpr auto index_ifft = 1;
      array<GrainType, n_non_multiplicative_grains> grain_types{{ GrainType::FFT, GrainType::IFFT }};
      array<float, n_non_multiplicative_grains> times;
      int index = 0;
      for(auto g : grain_types)
      {
          times[index] = min_(measure_thread_cpu_n(n_atoms_repeat_warmup,
                                                   n_atoms_repeat,
                                                   [&test, g] () { test.prepare(g);},
                                                   [&test   ] () { test.run();}));
        ++index;
      }

      struct PhasedCost : public Cost {
        GrainsCosts grains_costs;
          
          void logSubReport(std::ostream & os) const override {
              os << "grain fft  : " << grains_costs.fft << std::endl;
              os << "grain ifft : " << grains_costs.ifft << std::endl;
              os << "grain mult : " << grains_costs.mult << std::endl;
          }
      };

      struct CostEvaluator {
        array<float, n_non_multiplicative_grains> fft_times;
        int n_scales;
        int nAudioCbFrames;
        int n_channels;
        bool constraint;

        void evaluate(float multiplication_grain_time, int n_multiplicative_grains, int grain_period,
                      PhasedCost & result) const {
          result.setPhase(0);

          // factor 2 because the unit is half grain in 'grains_costs'.
          auto max_n_halfgrains_per_cb = 2 * nAudioCbFrames / grain_period;
          if(max_n_halfgrains_per_cb * grain_period != 2 * nAudioCbFrames) {
            ++max_n_halfgrains_per_cb; // worst case, not avg
          }

          int const n_half_grains = pow2(n_scales-1) * 2*(n_multiplicative_grains + n_non_multiplicative_grains);

          cyclic<float> grains_costs(n_half_grains);

          {
            // The order is:
            // fft, m1, m2, ... mn, ifft
            std::vector<float> costs;
            costs.push_back(fft_times[index_fft]);
            for(int i=0; i<n_multiplicative_grains; ++i) {
              costs.push_back(multiplication_grain_time);
            }
            costs.push_back(fft_times[index_ifft]);

            std::vector<int> scale_offset; // position of first computation in grains_costs
            for(int i=0; i<n_scales; ++i) {
              // ensures that for any number of scale, we will have
              // at most twice the density of computations for single scale
              scale_offset.push_back(pow2(i)-1);
           }
            for(int i=0; i<n_half_grains; ++i) {
              float cost = 0.f;
              for(int s=0; s<n_scales; ++s) {
                int ii = i-scale_offset[s];
                if(ii < 0) {
                  continue;
                }
                int div = ii / pow2(1+s);

                if(ii == div * pow2(1+s)) {
                  cost += costs[div % costs.size()];
                }
              }
              *(grains_costs.begin() + i) = cost;
            }
          }

          // TODO [early coefficients cost] we should have a sample-unit cyclic, and put one grain cost
          // every period

          float cost = computeMaxSlidingSum(grains_costs,
                                            max_n_halfgrains_per_cb);
          if(constraint) {
            assert(n_channels >= 2);
            int n_samples_between_grains = n_scales <= 1 ? grain_period : (grain_period/2);

            auto n_min_empty_cb_between_consecutive_grains = -1 + n_samples_between_grains / nAudioCbFrames;
            if(n_min_empty_cb_between_consecutive_grains >= n_channels - 1) {
              // easy case : there is enough room between grains to evenly distribute all channels
              result.setPhase(grain_period / n_channels);
            }
            else {
              // harder case: we need to go more in detail, and find the phase that minimizes
              // the worst callback cost.

              cost *= n_channels;
              // now cost is the 'phase == 0' cost
              result.setPhase(0);

              cyclic<float> phased_grains_costs(grains_costs.size());
              for(int phase = 1; phase < grains_costs.size(); ++phase) {
                compute_phased_sum(grains_costs,
                                   phase,
                                   n_channels,
                                   phased_grains_costs);

                auto phased_cost = computeMaxSlidingSum(phased_grains_costs,
                                                        max_n_halfgrains_per_cb);
                if(phased_cost < cost) {
                  cost = phased_cost;
                  result.setPhase(phase); // unit is "half of grain period"
                }
              }

              // convert phase units from "half of grain period" to "frames"
              result.setPhase(result.getPhase().value() * 2 * grain_period);
            }
            // cost is now the cost of each channel
            // but cost should be per sample, not per frame, so
            // we divide by the number of channels
            cost /= static_cast<float>(n_channels);
          }

          cost /= nAudioCbFrames;
          // cost == 'worst computation time over one callback, averaged per sample'

          result.grains_costs.fft  = fft_times[index_fft];
          result.grains_costs.ifft = fft_times[index_ifft];
          result.grains_costs.mult = multiplication_grain_time;

          result.setCost(cost);
        }
      } cost_evaluator{times, n_scales, n_frames, n_channels, constraint};

      RangedGradientDescent<PhasedCost> rgd([ &cost_evaluator, &test ](int multiplication_group_size, auto & cost) {
        // compute multiplication time for the group

        test.setMultiplicationGroupLength(multiplication_group_size);

        if(!test.isValid()) {
          return ParamState::OutOfRange;
        }
          auto multiplication_grain_time = min_(measure_thread_cpu_n(n_atoms_repeat_warmup,
                                                                     n_atoms_repeat,
                                                                     [&test](){ test.prepare(GrainType::MultiplicationGroup); },
                                                                     [&test](){ test.run(); }));

        // TODO [early coefficients cost] we should pass a vector of additional costs (where we take into account the early coefficient handling)

        cost_evaluator.evaluate(multiplication_grain_time,
                                test.countMultiplicativeGrains(),
                                test.getGranularMinPeriod(),
                                cost);

        //cout
        //<< "mult time for group size '" << multiplication_group_size << "' : " << multiplication_grain_time
        //<< " cost : '" << cost << "'" << endl;

        return ParamState::Ok;
      });

      range<int> const multiplication_group_length {
        test.getLowestValidMultiplicationsGroupSize(),
        test.getHighestValidMultiplicationsGroupSize()
      };

      PhasedCost phased_cost;
      val.multiplication_group_size = rgd.findLocalMinimum(n_iterations, multiplication_group_length, phased_cost);
      val.setCost(phased_cost.getCost());
      val.setGrainsCosts(phased_cost.grains_costs);
      val.setPhase(phased_cost.getPhase().value_or(0));
      val.partition_size = partition_size;

      constexpr auto debug = false;
      if(debug) {
        rgd.plot(true, os);
        rgd.make_exhaustive(multiplication_group_length, os);
        rgd.plot(true, os);
      }

      //            os
      //            << "optimal group size : " << val.multiplication_group_size
      //            << " cost : '" << cost << "'" << endl;

      return ParamState::Ok;
    });

    auto res = gradient_descent.findLocalMinimum(n_iterations,
                                                 min_lg2_partitionsz,
                                                 min_val);
    if(res) {
      assert(min_val);
      assert(min_val->partition_size == pow2(*res));
    }
  }

  template<typename NonAtomicConvolution, typename SetupParam = typename NonAtomicConvolution::SetupParam, typename GradientDescent = GradientDescent<SetupParam>>
  void get_optimal_partition_size_for_nonatomic_convolution(GradientDescent & gd,
                                                           int n_channels,
                                                           int n_scales,
                                                           bool with_spread,
                                                           int n_audiocb_frames,
                                                           std::function<Optional<int>(int)> n_coeffs_for_partition_sz,
                                                            int min_lg2_partition_sz,
                                                           Optional<SetupParam> & value,
                                                            std::ostream & os)
  {
    constexpr auto n_iterations = 1;
    find_optimal_partition_size_for_nonatomic_convolution<NonAtomicConvolution>(gd,
                                                                                n_iterations,
                                                                                n_channels,
                                                                                n_scales,
                                                                                n_audiocb_frames,
                                                                                n_coeffs_for_partition_sz,
                                                                                min_lg2_partition_sz,
                                                                                with_spread,
                                                                                value,
                                                                                os);
  }

  template<typename T>
  struct PartitionAlgo< FinegrainedPartitionnedFFTConvolution<T> > {
    using NonAtomicConvolution = FinegrainedPartitionnedFFTConvolution<T>;
    using SetupParam = typename NonAtomicConvolution::SetupParam;
    using PS = PartitionningSpec<SetupParam>;
    using PSpecs = PartitionningSpecs<SetupParam>;

    static PSpecs run(int n_channels,
                      int n_scales,
                      int n_audio_frames_per_cb,
                      std::function<Optional<int>(int)> n_coeffs_for_partition_sz,
                      int min_lg2_partition_sz,
                      std::ostream & os) {
      assert(n_channels > 0);
      PSpecs res;
      {
        auto & spec = res.without_spread;
        get_optimal_partition_size_for_nonatomic_convolution<NonAtomicConvolution>(spec.gd,
                                                                                   n_channels,
                                                                                   n_scales,
                                                                                   false,
                                                                                   n_audio_frames_per_cb,
                                                                                   n_coeffs_for_partition_sz,
                                                                                   min_lg2_partition_sz,
                                                                                   spec.cost,
                                                                                   os);
      }

      if(n_channels > 1) {
        auto & spec = res.with_spread;
        get_optimal_partition_size_for_nonatomic_convolution<NonAtomicConvolution>(spec.gd,
                                                                                   n_channels,
                                                                                   n_scales,
                                                                                   true,
                                                                                   n_audio_frames_per_cb,
                                                                                   n_coeffs_for_partition_sz,
                                                                                   min_lg2_partition_sz,
                                                                                   spec.cost,
                                                                                   os);
      }

      return std::move(res);
    }
  };


static std::ostream& operator<<(std::ostream& os, const GrainsCosts& g)
{
  using namespace std;
  os << "grain 'fft'  : " << g.fft << endl;
  os << "grain 'ifft' : " << g.ifft << endl;
  os << "grain 'mult' : " << g.mult << endl;
  return os;
}
  static inline std::ostream& operator<<(std::ostream& os, const FinegrainedSetupParam& p)
  {
    using namespace std;
    os
    << p.grains_costs
    << "multiplication group size : " << p.multiplication_group_size << endl;
    return os;
  }
}
