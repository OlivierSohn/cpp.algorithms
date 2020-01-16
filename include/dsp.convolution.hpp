
namespace imajuscule
{
    template <typename Parent>
    struct FFTConvolutionIntermediateSimulation : public Parent {
        using FPT = typename Parent::FPT;
        using FFTTag = typename Parent::FFTTag;
        using FFTAlgo = typename fft::Algo_<FFTTag, FPT>;
        using Parent::doSetCoefficientsCount;
        static constexpr bool has_subsampling = false;

        static constexpr auto cost_copy = fft::RealSignalCosts<FFTTag, FPT>::cost_copy;
        static constexpr auto cost_add_scalar_multiply = fft::RealSignalCosts<FFTTag, FPT>::cost_add_scalar_multiply;
        static constexpr auto cost_fft_inverse = fft::AlgoCosts<FFTTag, FPT>::cost_fft_inverse;
        
        using Parent::cost_compute_convolution;
        using Parent::get_fft_length;
        using Parent::getBlockSize;
        
        void setCoefficientsCount(int64_t szCoeffs) {
            doSetCoefficientsCount(szCoeffs);
        }

        FFTConvolutionIntermediateSimulation()
        : minorCost(costWriteNConsecutive<FPT*>(1) +
                    costReadNConsecutive<FPT>(1))
        {}
        
        double simuMinorStep() {
            return minorCost;
        }

        double simuMajorStep() {
            if(unlikely(!majorCost)) {
                majorCost =
                cost_compute_convolution() + // assuming that cost_compute_convolution is CONSTANT;
                cost_fft_inverse(get_fft_length()) +
                cost_add_scalar_multiply(getBlockSize()) +
                cost_copy(getBlockSize()) +
                costReadNConsecutive<FPT>(1);
            }
            return *majorCost;
        }
    private:
        double minorCost;
        std::optional<double> majorCost;
    };

    /*
     * cf. https://en.wikipedia.org/wiki/Overlap%E2%80%93add_method
     */
    template <typename Parent>
    struct FFTConvolutionIntermediate : public Parent {
        using Simulation = FFTConvolutionIntermediateSimulation<typename Parent::Simulation>;
        
        using T = typename Parent::FPT;
        using FPT = T;
        using Tag = typename Parent::FFTTag;
        static constexpr int nCoefficientsFadeIn = 0;
        static constexpr bool has_subsampling = false;
        static constexpr bool step_can_error = false;

      using SetupParam = typename Parent::SetupParam;

        static constexpr auto copy = fft::RealSignal_<Tag, FPT>::copy;
        static constexpr auto get_signal = fft::RealSignal_<Tag, FPT>::get_signal;
        static constexpr auto zero_signal = fft::RealSignal_<Tag, FPT>::zero;
        static constexpr auto add_scalar_multiply = fft::RealSignal_<Tag, FPT>::add_scalar_multiply;
        using RealSignal = typename fft::RealSignal_<Tag, FPT>::type;
        using Algo = typename fft::Algo_<Tag, FPT>;
        using Contexts = fft::Contexts_<Tag, FPT>;

        using Parent::get_fft_length;
      using Parent::getBlockSize;
      using Parent::getEpsilon;
      using Parent::doLogComputeState;
      using Parent::doSetCoefficients;
        using Parent::doFlushToSilence;
        using Parent::compute_convolution;

        FFTConvolutionIntermediate() {
            it = y.end();
        }
        void logComputeState(std::ostream & os) const {
            os << "FFTConv [" << std::distance(y.begin(), it) << "/" << getBlockSize() << "]" << std::endl;
            doLogComputeState(os);
        }

        void setCoefficients2(a64::vector<T> coeffs_) {
            auto const N = coeffs_.size();
            auto const fft_length = get_fft_length(N);
            fft.setContext(Contexts::getInstance().getBySize(fft_length));

            result.clear();
            result.resize(fft_length);

            y.clear();
            y.resize(fft_length);
            it = y.begin();

            doSetCoefficients(fft, std::move(coeffs_));
        }
      
      void reset() {
        y.clear();
        it = y.end();
        result.clear();
      }
      void flushToSilence() {
        Parent::doFlushToSilence();
        zero_signal(y);
        zero_signal(result);
        it = y.begin();
      }
      
      bool isZero() const {
        return result.empty();
      }

        // minor step
      FPT doStep() {
        ++it;
        assert(it < y.begin() + getBlockSize());
        return get_signal(*it);
      }

        // major step
      // the input vector is expected to be padded.
      FPT doStep(typename RealSignal::const_iterator xBegin)
      {
        auto const N = getBlockSize();
        auto const & frequencies = compute_convolution(fft, xBegin);

        fft.inverse(frequencies, result, get_fft_length());

        auto factor = 1 / (Algo::scale * Algo::scale * static_cast<T>(get_fft_length()));

        auto it_res = result.begin();
        auto it_y = y.begin();
        auto it_y_prev = it_y + N;

        // y = mix first part of result with second part of previous result
        //
        // 'first part of y' = factor * ('second part of y' + 'first part of result')
        add_scalar_multiply(it_y, /* = */
                            /* ( */ it_res, /* + */ it_y_prev /* ) */, /* x */ factor,
                            N);

        // store second part of result for later
        //
        // 'second part of y' = 'second part of result'
        copy(it_y   + N,
             it_res + N,
             N);

        assert(it == y.begin() + N-1); // make sure 'rythm is good', i.e we exhausted the first half of the y vector
        it = y.begin();
        return get_signal(*it);
      }

    private:
        Algo fft;
        RealSignal y, result;
        typename decltype(y)::const_iterator it;
    };

  template <typename Parent>
  struct FFTConvolutionBase : public Parent {
    using T = typename Parent::FPT;
    using FPT = T;
    using Tag = typename Parent::FFTTag;
      using Parent::flushToSilence;
    static constexpr int nCoefficientsFadeIn = 0;

    using SetupParam = typename Parent::SetupParam;
    
    using RealSignal = typename fft::RealSignal_<Tag, FPT>::type;
      static constexpr auto zero_signal = fft::RealSignal_<Tag, FPT>::zero;

    using Parent::get_fft_length;
    using Parent::getBlockSize;
    using Parent::getEpsilon;
    using Parent::doStep;
    using Parent::setCoefficients2;

    void setCoefficients(a64::vector<T> coeffs_) {
        progress = 0;
      x.clear();
      auto const fft_length = get_fft_length(coeffs_.size());
      x.resize(fft_length); // padding
      setCoefficients2(std::move(coeffs_));
    }

      void flushToSilence() {
          Parent::flushToSilence();
          progress = 0;
      }

    bool willComputeNextStep() const {
      return progress == getBlockSize()-1;
    }

    T step(T val) {
      x[progress] = typename RealSignal::value_type(val);
      ++progress;

      if(progress == getBlockSize()) {
        // x is already padded
        auto res = doStep(x.begin());
        progress = 0;
        return res;
      }
      else {
        return doStep();
      }
    }
      
      template<typename FPT2>
      void stepAddVectorized(FPT2 const * const input_buffer,
                             FPT2 * output_buffer,
                             int nSamples)
      {
          for(int i=0; i<nSamples; ++i) {
              output_buffer[i] += step(input_buffer[i]);
          }
      }

  private:
    RealSignal x;
    int progress = 0;
  };

struct FFTConvolutionCRTPSetupParam : public Cost
{
    FFTConvolutionCRTPSetupParam(int blockSize)
    : blockSize(blockSize)
    {}
    
    void logSubReport(std::ostream & os) const override {
        os << "FFTConvolutionCRTPSetupParam block size " << blockSize << std::endl;
    }
    
    int getImpliedLatency() const {
        return blockSize-1;
    }

    int blockSize;
};

template <typename T, typename Tag>
struct FFTConvolutionCRTPSimulation {
    using SetupParam = FFTConvolutionCRTPSetupParam;
    using FPT = T;
    using FFTTag = Tag;
    using FFTAlgo = typename fft::Algo_<Tag, FPT>;
    
    static constexpr auto cost_mult_assign = fft::RealFBinsCosts<Tag, FPT>::cost_mult_assign;
    static constexpr auto cost_fft_forward = fft::AlgoCosts<Tag, FPT>::cost_fft_forward;

    void doSetCoefficientsCount(int64_t szCoeffs) {
        fft_length = get_fft_length(szCoeffs);
        N = fft_length/2;
    }

    double cost_compute_convolution() {
        return
        cost_fft_forward(get_fft_length()) +
        cost_mult_assign(get_fft_length());
    }
    
    static auto get_fft_length(int n) {
        auto N_nonzero_y = 2 * n;
        return ceil_power_of_two(N_nonzero_y);
    }

    auto get_fft_length() const {
        return fft_length;
    }
    
    auto getBlockSize() const { return N; }

    int getLatency() const { return N-1; }

    int N, fft_length = 0;
};

    template <typename T, typename Tag>
    struct FFTConvolutionCRTP {
        using Simulation = FFTConvolutionCRTPSimulation<T, Tag>;
        
        using FPT = T;
        using FFTTag = Tag;

        using RealSignal = typename fft::RealSignal_<Tag, FPT>::type;
        static constexpr auto makeRealSignal = fft::RealSignal_<Tag, FPT>::make;

        using CplxFreqs = typename fft::RealFBins_<Tag, FPT>::type;
        static constexpr auto mult_assign = fft::RealFBins_<Tag, FPT>::mult_assign;

        using Algo = typename fft::Algo_<Tag, FPT>;

        using SetupParam = FFTConvolutionCRTPSetupParam;
      
      void setup(SetupParam const & p) {
          N = p.blockSize;
      }
        
        void doLogComputeState(std::ostream & os) const {
            os << "1 block of size " << N << std::endl;
        }
    private:
        int N;

        CplxFreqs fft_of_h, fft_of_x;
    public:

        auto getBlockSize() const { return N; }
        auto getLatency() const { return N-1; }
        auto getGranularMinPeriod() const { return getBlockSize(); }
        bool isValid() const { return true; }

        auto countPartitions() const { return 1; }

        bool empty() const { return fft_of_h.empty(); }

        static auto get_fft_length(int n) {
            auto N_nonzero_y = 2 * n;
            return ceil_power_of_two(N_nonzero_y);
        }

        auto get_fft_length() const {
            return fft_of_h.size();
        }

        void doSetCoefficients(Algo const & fft, a64::vector<T> coeffs_) {
            auto fft_length = get_fft_length(coeffs_.size());
            Assert(N == fft_length/2);

            fft_of_h.resize(fft_length);
            fft_of_x.resize(fft_length);

            // pad impulse response with 0

            coeffs_.resize(fft_length, {});

            // compute fft of padded impulse response
            auto coeffs = makeRealSignal(std::move(coeffs_));
            fft.forward(coeffs.begin(), fft_of_h, fft_length);
        }


      double getEpsilon() const {
        return fft::getFFTEpsilon<FPT>(get_fft_length()) + 2 * std::numeric_limits<FPT>::epsilon();
      }

    protected:

      auto const & compute_convolution(Algo const & fft, typename RealSignal::const_iterator x) {
            // do fft of x

            fft.forward(x, fft_of_x, get_fft_length());

            mult_assign(fft_of_x, fft_of_h);

            return fft_of_x;
        }
        
        void doFlushToSilence() {
            // nothing to do : no member contains state related to a past signal, except
            // for fft_of_x but it will be overwritten next time we use it in compute_convolution
        }
    };

    /*
     * Partitionned convolution, cf. http://www.ericbattenberg.com/school/partconvDAFx2011.pdf
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
     * Notes : Convolution in "time" space is the same as multiplication in "frequency" space.
     */

struct PartitionnedFFTConvolutionCRTPSetupParam {
    PartitionnedFFTConvolutionCRTPSetupParam(int partition_size)
    :partition_size(partition_size)
    {}
    
    int getImpliedLatency() const {
        return partition_size-1;
    }
    int partition_size;
    
    bool isValid() const {
        return partition_size > 0;
    }
    
    void logSubReport(std::ostream & os) const {
        os << "partition_sz : " << partition_size << std::endl;
    }
};

template <typename T, typename Tag>
struct PartitionnedFFTConvolutionCRTPSimulation {
    using SetupParam = PartitionnedFFTConvolutionCRTPSetupParam;
    using FPT = T;
    using FFTTag = Tag;

    static constexpr auto cost_multiply     = fft::RealFBinsCosts<Tag, FPT>::cost_multiply;
    static constexpr auto cost_multiply_add = fft::RealFBinsCosts<Tag, FPT>::cost_multiply_add;
    static constexpr auto cost_fft_forward = fft::AlgoCosts<Tag, FPT>::cost_fft_forward;

    void setup(SetupParam const & p) {
      partition_size = p.partition_size;
      assert(partition_size > 0);
      assert(is_power_of_two(partition_size));
    }

    void doSetCoefficientsCount(int64_t count) {
        n_partitions = [&count, partition_size = this->partition_size](){
            auto const N = count;
            auto n_partitions = N/partition_size;
            if(n_partitions * partition_size != N) {
                // one partition is partial...
                assert(n_partitions * partition_size < N);
                ++n_partitions;
                // ... pad it with zeros
                count = n_partitions * partition_size;
            }
            return n_partitions;
        }();
    }
    
    double cost_compute_convolution() const
    {
        double cost = cost_fft_forward(get_fft_length());
        if(n_partitions) {
            cost += cost_multiply(get_fft_length());
            cost += (n_partitions-1) * cost_multiply_add(get_fft_length());
        }
        return cost;
    }
    
    auto get_fft_length() const { assert(partition_size > 0); return 2 * partition_size; }
    auto getBlockSize() const { assert(partition_size > 0); return partition_size; }
    auto getLatency() const { assert(partition_size > 0); return partition_size-1; }
private:
    int partition_size = -1;
    int n_partitions = 0;
};

    template <typename T, typename Tag>
    struct PartitionnedFFTConvolutionCRTP {
        using Simulation = PartitionnedFFTConvolutionCRTPSimulation<T, Tag>;
        using SetupParam = PartitionnedFFTConvolutionCRTPSetupParam;
        using FPT = T;
        using FFTTag = Tag;

        using RealSignal = typename fft::RealSignal_<Tag, FPT>::type;
        using Signal_value_type = typename fft::RealSignal_<Tag, FPT>::value_type;

        using CplxFreqs = typename fft::RealFBins_<Tag, FPT>::type;
        static constexpr auto zero = fft::RealFBins_<Tag, FPT>::zero;
        static constexpr auto multiply = fft::RealFBins_<Tag, FPT>::multiply;
        static constexpr auto multiply_add = fft::RealFBins_<Tag, FPT>::multiply_add;

        using Algo = typename fft::Algo_<Tag, FPT>;
        
        void doLogComputeState(std::ostream & os) const {
            os << countPartitions() << " partitions " << std::endl;
        }

        auto get_fft_length() const { assert(partition_size > 0); return 2 * partition_size; }
        auto get_fft_length(int) const { return get_fft_length(); }

        bool empty() const { return ffts_of_partitionned_h.empty(); }

        auto getBlockSize() const { return partition_size; }
        auto getGranularMinPeriod() const { return getBlockSize(); }
        auto getLatency() const { return partition_size-1; }

        bool isValid() const { return true; }

        auto countPartitions() const { return ffts_of_partitionned_h.size(); }

      

      void setup(SetupParam const & p) {
        partition_size = p.partition_size;
        assert(partition_size > 0);
        assert(is_power_of_two(partition_size));
      }

        void doSetCoefficients(Algo const & fft, a64::vector<T> coeffs_) {

            auto const n_partitions = [&coeffs_, partition_size = this->partition_size](){
                auto const N = coeffs_.size();
                auto n_partitions = N/partition_size;
                if(n_partitions * partition_size != N) {
                    // one partition is partial...
                    assert(n_partitions * partition_size < N);
                    ++n_partitions;
                    // ... pad it with zeros
                    coeffs_.resize(n_partitions * partition_size, {});
                }
                return n_partitions;
            }();

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

      double getEpsilon() const {
        return countPartitions() * (fft::getFFTEpsilon<FPT>(get_fft_length()) + 2 * std::numeric_limits<FPT>::epsilon());
      }
    protected:

        void doFlushToSilence() {
            for(auto & fft_of_delayed_x : ffts_of_delayed_x) {
                zero(fft_of_delayed_x);
            }
            ffts_of_delayed_x.setByIndex(0);
        }

        auto const & compute_convolution(Algo const & fft, typename RealSignal::const_iterator & xBegin)
        {
            // do fft of x
            {
                auto & oldest_fft_of_delayed_x = *ffts_of_delayed_x.cycleEnd();
                auto const fft_length = get_fft_length();
                assert(fft_length == oldest_fft_of_delayed_x.size());
                fft.forward(xBegin, oldest_fft_of_delayed_x, fft_length);
                ffts_of_delayed_x.advance();
            }

            int index = 0;

            ffts_of_delayed_x.for_each_bkwd( [this, &index] (auto const & fft_of_delayed_x) {
                assert(index < ffts_of_partitionned_h.size());
                auto const & fft_of_partitionned_h = ffts_of_partitionned_h[index];
                
                if(index == 0) {
                    multiply(work /* = */, fft_of_delayed_x, /* * */ fft_of_partitionned_h);
                }
                else {
                    multiply_add(work /* += */, fft_of_delayed_x, /* * */ fft_of_partitionned_h);
                }
                ++ index;
            });

            return work;
        }

    private:
        int partition_size = -1;
        cyclic<CplxFreqs> ffts_of_delayed_x;
        std::vector<CplxFreqs> ffts_of_partitionned_h;

        CplxFreqs work;
    };

    /*

     Notations for complexity:

     H : length of impulse response

     A : number of frames computed during one audio callback

     PART_N : Number of partitions
     PART_L : Length of one partition
     ( H == PART_N * PART_L )

     The computations are based on the fact that an fft of an S-long signal costs S*lg(S)
     */

    /*
     * runtime complexity:
     *
     *   every H frames  ............... O( H * lg(H) )
     *
     *   every frame  .................. O( lg(H) )       [amortized]
     *
     *   when 'H < A':
     *     worst audio callback call ... O( A * lg(H) )    (= A/H * O( H * lg(H) ))
     *
     *   when 'A < H ':
     *     worst audio callback call ... O( H * lg(H) )
     *
     * optimization : H and A are fixed so we cannot optimize this algorithm
     */
  template <typename T, typename FFTTag>
  using FFTConvolutionCore = FFTConvolutionIntermediate < FFTConvolutionCRTP<T, FFTTag> >;
  template <typename T, typename FFTTag>
  using FFTConvolution = FFTConvolutionBase< FFTConvolutionCore<T, FFTTag> >;

    /*
     * runtime complexity:
     *
     *   every PART_L frames                :   O( PART_L * (PART_N + lg(PART_L) ) )
     *
     *   every frame                        :   O( PART_N + lg(PART_L) )      [amortized]
     *
     *   with 'PART_L < A':
     *     worst audio callback call ... O(     A  * (part_N + lg(PART_L)) )       (= A/PART_L *  O( PART_L * (PART_N + lg(PART_L) ) )
     *
     *   with 'A < PART_L '
     *     worst audio callback call ... O( PART_L * (PART_N + lg(PART_L) ) )
     *                                 = O( PART_L * (H/PART_L + lg(PART_L) ) )
     *                                 = O( H + PART_L * lg(PART_L) ) )
     *
     * optimization : PART_L is not fixed so we can optimize this algorithm by trying different powers of 2
     *                also we can optimize more globally, taking into account that we have one reverb per channel:
     *                when N_Channel * A < PART_L we can distribute the computes over the different callbacks calls, provided the "phases"
     *                of the different algorithms are well-spaced.
     *
     */
    template <typename T, typename FFTTag>
    using PartitionnedFFTConvolution = FFTConvolutionBase< FFTConvolutionIntermediate < PartitionnedFFTConvolutionCRTP<T, FFTTag> > >;
}
