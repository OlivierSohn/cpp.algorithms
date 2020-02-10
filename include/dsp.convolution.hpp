
namespace imajuscule {

template <typename Parent>
struct FFTConvolutionIntermediateSimulation : public Parent {
    using FPT = typename Parent::FPT;
    using FFTTag = typename Parent::FFTTag;
    using Parent::doSetCoefficientsCount;
    static constexpr bool has_subsampling = false;
    
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
    
    double simuMajorStep(XFFTsCostsFactors const & xFftCostFactors) {
        if(unlikely(!majorCostExceptXForwardFft)) {
            majorCostExceptXForwardFft =
            
            // assuming that cost_compute_convolution is CONSTANT
            cost_compute_convolution() +
            
            fft::AlgoCosts<FFTTag, FPT>::cost_fft_inverse(get_fft_length()) +
            
            // if get_fft_length() has the same order of magnitude than the size of y:
            // - many 'add_assign' will be replaced by 'copy'
            // - many 'add_assign' will be done in 2 separate passes
            //
            // To take these effects fully into account we would need to know
            // - the size of y
            // - the phase at which we will write in y
            // - the rythm of other convolution parts (and other channels) writing in the same y
            //
            // "This is complicated to do" and we tend to think that the magnitude of the error
            // will be small compared to the overall sum, hence we simulate using a single 'add_assign':
            
            fft::RealSignalCosts<FFTTag, FPT>::cost_add_assign(get_fft_length()) +
            
            costReadNConsecutive<FPT>(1);
        }
        if(unlikely(!costXForwardFft)) {
            costXForwardFft = fft::AlgoCosts<FFTTag, FPT>::cost_fft_forward(get_fft_length());
        }
        //LG(INFO, "*:%f fft: %f", *majorCostExceptXForwardFft, *costXForwardFft);
        return
        *majorCostExceptXForwardFft +
        *costXForwardFft * xFftCostFactors.get(get_fft_length());
    }
private:
    double minorCost;
    std::optional<double> majorCostExceptXForwardFft;
    std::optional<double> costXForwardFft;
};


/*
 * cf. https://en.wikipedia.org/wiki/Overlap%E2%80%93add_method
 */
template <typename Parent>
struct FFTConvolutionIntermediate : public Parent {
    
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
    
    template<typename X>
    using Allocator = typename Parent::template Allocator<X>;
    
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
        
        fft.inverse(frequencies.data(), result.data(), get_fft_length());
        
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
        copy(&y[N],
             &result[N],
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
    static constexpr bool has_subsampling = false;
    static constexpr int nCoefficientsFadeIn = 0;

    FFTConvolutionCRTPSetupParam(int blockSize)
    : blockSize(blockSize)
    {}
    
    void logSubReport(std::ostream & os) const override {
        os << "FFTConvolutionCRTPSetupParam block size " << blockSize << std::endl;
    }
    
    bool handlesCoefficients() const {
        return blockSize > 0;
    }
    
    template<Overlap Mode>
    MinSizeRequirement getMinSizeRequirement() const
    {
        int const fft_length = 2 * blockSize;
        
        int const y_size = [fft_length](){
            if constexpr(Mode == Overlap::Add) {
                return fft_length;
            }
            else {
                return fft_length/2;
            }
        }();
        
        return {
            0, // x block size
            y_size, // y block size
            {
                {fft_length, 1}
            },
            2*fft_length // work size
        };
    }
    
    Latency getImpliedLatency() const {
        Assert(handlesCoefficients());
        return Latency(blockSize-1);
    }

    template<typename F>
    void forEachUsingSameContext(F f) const {
        f(*this);
    }

    int blockSize;
};

template <typename T, typename Tag>
struct FFTConvolutionCRTPSimulation {
    using SetupParam = FFTConvolutionCRTPSetupParam;
    using FPT = T;
    using FFTTag = Tag;
    
    void doSetCoefficientsCount(int64_t szCoeffs) {
        fft_length = get_fft_length(szCoeffs);
        N = fft_length/2;
    }

    double cost_compute_convolution() {
        return fft::RealFBinsCosts<Tag, FPT>::cost_mult_assign(N);
    }
    
    static auto get_fft_length(int n) {
        auto N_nonzero_y = 2 * n;
        return ceil_power_of_two(N_nonzero_y);
    }

    auto get_fft_length() const {
        return fft_length;
    }
    
    auto getBlockSize() const { return N; }

    bool handlesCoefficients() const {
        return N > 0;
    }
    Latency getLatency() const {
        Assert(handlesCoefficients());
        return Latency(N-1);
    }

    int N, fft_length = 0;
};

template<typename T, typename FFTTag>
struct Simulation_<FFTConvolutionCRTPSetupParam, T, FFTTag> {
    using type = FFTConvolutionIntermediateSimulation<FFTConvolutionCRTPSimulation<T, FFTTag>>;
};

template <typename T, template<typename> typename Alloc, typename Tag>
struct FFTConvolutionCRTP {
    template<typename TT>
    using Allocator = Alloc<TT>;
    
    using FPT = T;
    using FFTTag = Tag;
    
    using RealSignal = typename fft::RealSignal_<Tag, FPT>::type;
    static constexpr auto makeRealSignal = fft::RealSignal_<Tag, FPT>::make;
    
    using CplxFreqs = typename fft::RealFBins_<Tag, FPT, Allocator>::type;
    static constexpr auto mult_assign = fft::RealFBins_<Tag, FPT, Allocator>::mult_assign;
    
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
    bool handlesCoefficients() const {
        return getBlockSize() > 0;
    }
    Latency getLatency() const {
        Assert(handlesCoefficients());
        return Latency(N-1);
    }
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
        fft.forward(coeffs.begin(), fft_of_h.data(), fft_length);
    }
    
    
    double getEpsilon() const {
        return fft::getFFTEpsilon<FPT>(get_fft_length()) + 2 * std::numeric_limits<FPT>::epsilon();
    }
    
protected:
    
    auto const & compute_convolution(Algo const & fft, typename RealSignal::const_iterator x) {
        // do fft of x
        
        fft.forward(x, fft_of_x.data(), get_fft_length());
        
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

struct PartitionnedFFTConvolutionCRTPSetupParam : public Cost {
    static constexpr bool has_subsampling = false;
    static constexpr int nCoefficientsFadeIn = 0;

    PartitionnedFFTConvolutionCRTPSetupParam(int partition_size,
                                             int partition_count)
    : partition_size(partition_size)
    , partition_count(partition_count)
    {}
    
    void adjustWork(int const nTargetCoeffs) {
        partition_count = countPartitions(nTargetCoeffs,
                                          partition_size);
        if(!handlesCoefficients()) {
            setCost(0.);
        }
    }
    
    template<Overlap Mode>
    MinSizeRequirement getMinSizeRequirement() const
    {
        int const fft_length = 2 * partition_size;
        
        int const y_size = [fft_length](){
            if constexpr(Mode == Overlap::Add) {
                return fft_length;
            }
            else {
                return fft_length/2;
            }
        }();
        
        return {
            0, // x block size
            y_size, // y block size
            {
                {fft_length, partition_count}
            },
            2*fft_length // work size
        };
    }
    
    bool handlesCoefficients() const {
        return countMaxHandledCoeffs() > 0;
    }
    int countMaxHandledCoeffs() const {
        return partition_count * partition_size;
    }
    Latency getImpliedLatency() const {
        Assert(handlesCoefficients());
        return Latency(partition_size-1);
    }
    
    int partition_size;
    int partition_count;
    
    bool isValid() const {
        return partition_size > 0;
    }
    
    template<typename F>
    void forEachUsingSameContext(F f) const {
        f(*this);
    }

    void logSubReport(std::ostream & os) const {
        os << partition_count << " partition(s) of size " << partition_size << std::endl;
    }
};

template <typename T, typename Tag>
struct PartitionnedFFTConvolutionCRTPSimulation {
    using SetupParam = PartitionnedFFTConvolutionCRTPSetupParam;
    using FPT = T;
    using FFTTag = Tag;

    void setup(SetupParam const & p) {
      partition_size = p.partition_size;
      assert(partition_size > 0);
      assert(is_power_of_two(partition_size));
    }

    void doSetCoefficientsCount(int64_t count) {
        n_partitions = countPartitions(count,
                                       partition_size);
    }
    
    double cost_compute_convolution() const
    {
        double cost = 0.;
        if(n_partitions) {
            cost += fft::RealFBinsCosts<Tag, FPT>::cost_multiply(partition_size);
            cost += (n_partitions-1) * fft::RealFBinsCosts<Tag, FPT>::cost_multiply_add(partition_size);
        }
        return cost;
    }
    
    auto get_fft_length() const { assert(partition_size > 0); return 2 * partition_size; }
    auto getBlockSize() const { assert(partition_size > 0); return partition_size; }
    
    bool handlesCoefficients() const {
        return n_partitions * partition_size > 0;
    }
    Latency getLatency() const {
        Assert(handlesCoefficients());
        return Latency(partition_size-1);
    }
private:
    int partition_size = -1;
    int n_partitions = 0;
};

template<typename T, typename FFTTag>
struct Simulation_<PartitionnedFFTConvolutionCRTPSetupParam, T, FFTTag> {
    using type = FFTConvolutionIntermediateSimulation<PartitionnedFFTConvolutionCRTPSimulation<T, FFTTag>>;
};


template <typename T, template<typename> typename Alloc, typename Tag>
struct PartitionnedFFTConvolutionCRTP {
    template<typename TT>
    using Allocator = Alloc<TT>;
    using SetupParam = PartitionnedFFTConvolutionCRTPSetupParam;
    using FPT = T;
    using FFTTag = Tag;
    
    using RealSignal = typename fft::RealSignal_<Tag, FPT>::type;
    using Signal_value_type = typename fft::RealSignal_<Tag, FPT>::value_type;
    
    using RealFBins = typename fft::RealFBins_<Tag, FPT, Allocator>;
    using CplxFreqs = typename RealFBins::type;
    static constexpr auto zero = RealFBins::zero;
    
    using Algo = typename fft::Algo_<Tag, FPT>;
    
    void doLogComputeState(std::ostream & os) const {
        os << countPartitions() << " partitions " << std::endl;
    }
    
    auto get_fft_length() const { assert(partition_size > 0); return 2 * partition_size; }
    auto get_fft_length(int) const { return get_fft_length(); }
    
    bool empty() const { return ffts_of_partitionned_h.empty(); }
    
    auto getBlockSize() const { return partition_size; }
    auto getGranularMinPeriod() const { return getBlockSize(); }
    bool handlesCoefficients() const {
        return getBlockSize() * countPartitions() > 0;
    }
    Latency getLatency() const {
        Assert(handlesCoefficients());
        return Latency(partition_size-1);
    }
    
    bool isValid() const { return true; }
    
    int countPartitions() const { return partition_count; }
    
    static int getAllocationSz_SetCoefficients(SetupParam const & p) {
        int const fft_length = 2 * p.partition_size;
        return fft_length * (2 * p.partition_count);
    }
    
    void setup(SetupParam const & p) {
        partition_size = p.partition_size;
        partition_count = p.partition_count;
        assert(partition_size > 0);
        assert(is_power_of_two(partition_size));
    }
    
    void doSetCoefficients(Algo const & fft, a64::vector<T> coeffs_) {
        auto const fft_length = get_fft_length();
        
        work.resize(fft_length);
        
        auto const n_partitions = imajuscule::countPartitions(coeffs_.size(),
                                                              partition_size);
        if(n_partitions != countPartitions()) {
            throw std::logic_error("inconsistent number of partitions");
        }
        // if one partition is partial, pad it with zeros
        coeffs_.resize(n_partitions * partition_size,
                       {});
        
        ffts_of_delayed_x.resize(n_partitions);
        ffts_of_partitionned_h.resize(n_partitions);
        
        for(auto & fft_of_partitionned_h : ffts_of_partitionned_h) {
            fft_of_partitionned_h.resize(fft_length);
        }
        for(auto & fft_of_delayed_x : ffts_of_delayed_x) {
            fft_of_delayed_x.resize(fft_length);
        }
        
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
                fft.forward(coeffs_slice.begin(), fft_of_partitionned_h.data(), fft_length);
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
        if(!ffts_of_delayed_x.empty()) {
            ffts_of_delayed_x.setByIndex(0);
        }
    }
    
    auto const & compute_convolution(Algo const & fft, typename RealSignal::const_iterator & xBegin)
    {
        // do fft of x
        {
            ffts_of_delayed_x.go_back(); // we write backwards so that we can traverse the container forwards for multiply_adds
            auto & oldest_fft_of_delayed_x = *ffts_of_delayed_x.cycleEnd();
            auto const fft_length = get_fft_length();
            assert(fft_length == oldest_fft_of_delayed_x.size());
            fft.forward(xBegin, oldest_fft_of_delayed_x.data(), fft_length);
        }
        
        int index = 0;
        
        ffts_of_delayed_x.for_each( [this, &index] (auto const & fft_of_delayed_x) {
            assert(index < ffts_of_partitionned_h.size());
            auto const & fft_of_partitionned_h = ffts_of_partitionned_h[index];
            
            if(index == 0) {
                RealFBins::multiply(    work.data() /*  = */, fft_of_delayed_x.data(), /* * */ fft_of_partitionned_h.data(),
                                    partition_size);
            }
            else {
                RealFBins::multiply_add(work.data() /* += */, fft_of_delayed_x.data(), /* * */ fft_of_partitionned_h.data(),
                                        partition_size);
            }
            ++ index;
        });
        
        return work;
    }
    
private:
    int partition_size = -1;
    int partition_count = 0;
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
template <typename T, template<typename> typename Allocator, typename FFTTag>
using FFTConvolutionCore = FFTConvolutionIntermediate < FFTConvolutionCRTP<T, Allocator, FFTTag> >;
template <typename T, template<typename> typename Allocator, typename FFTTag>
using FFTConvolution = FFTConvolutionBase< FFTConvolutionCore<T, Allocator, FFTTag> >;

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
template <typename T, template<typename> typename Allocator, typename FFTTag>
using PartitionnedFFTConvolution =
FFTConvolutionBase< FFTConvolutionIntermediate < PartitionnedFFTConvolutionCRTP<T, Allocator, FFTTag> > >;
}
