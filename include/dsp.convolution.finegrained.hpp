
namespace imajuscule {

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
 * An algorithm computes the optimal parameters (to minimize the worst case cost for a single callback):
 *  - The size of the partitions
 *  - The count of multiplications per multiplication grain
 *  - The "phasing" of simultaneous convolutions to best interleave
 *      high-cost grains.
 * based on :
 *  - The callback buffer size
 *  - The count of simultaneous convolutions happening in a callback
 */

enum class GrainType {
    FFT,
    IFFT,
    MultiplicationGroup
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
    static constexpr int nCoefficientsFadeIn = 0;
    static constexpr bool has_subsampling = false;

    explicit FinegrainedSetupParam() {}
    
    FinegrainedSetupParam(int partitionSz,
                          int partition_count,
                          int multiplication_group_size,
                          int phase)
    : Cost(phase)
    , multiplication_group_size(multiplication_group_size)
    , partition_size(partitionSz)
    , partition_count(partition_count)
    {}
    
    void setGrainsCosts(GrainsCosts gcosts) { grains_costs = gcosts; }
    
    int multiplication_group_size = 0;
    int partition_size = 0;
    int partition_count = 0;
    GrainsCosts grains_costs;
    
    void logSubReport(std::ostream & os) const override {
        os << "Finegrained, " << partition_count << " partitions of size " << partition_size << ", mult group size: " << multiplication_group_size << std::endl;
    }
    
    bool handlesCoefficients() const {
        return partition_size > 0;
    }
    
    void adjustWork(int targetNCoeffs) {
        partition_count = countPartitions(targetNCoeffs, partition_size);
        if(!handlesCoefficients()) {
            setCost(0.);
        }
        else {
            // we could reduce the cost...
        }
    }
    
    int get_fft_length() const {
        return 2 * partition_size;
    }

    int getGranularity() const {
        return partition_size / countGrains();
    }
    
    static constexpr int countNonMultiplicativeGrains() { return 2; }

    int countMultiplicativeGrains() const
    {
        return multiplication_group_size ?
        (1 + (partition_count - 1) / multiplication_group_size) :
        0;
    }
    
    int countGrains() const {
        return countMultiplicativeGrains() + countNonMultiplicativeGrains();
    }
    
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
        if(partition_count == 0) {
            return 0;
        }
        for(int i=1;; ++i) {
            if( (partition_count - 1)/i <= diff) {
                return i;
            }
        }
    }
    
    int getHighestValidMultiplicationsGroupSize() const {
        return partition_count;
    }
    
    MinSizeRequirement getMinSizeRequirement() const
    {
        int const blockProgressForIFFTGrain = (countMultiplicativeGrains()+1) * getGranularity();
        Assert(blockProgressForIFFTGrain <= partition_size);
        int const gap = partition_size - blockProgressForIFFTGrain;

        return {
            0, // x block size
            static_cast<int>(get_fft_length() + gap), // y block size
            {
                {get_fft_length(), partition_count}
            },
            get_fft_length() // work size
        };
    }
    
    static constexpr Latency getLatencyForPartitionSize(int sz) {
        return Latency(2 * sz - 1);
    }

    Latency getImpliedLatency() const {
        Assert(handlesCoefficients());
        return getLatencyForPartitionSize(partition_size);
    }
    
    template<typename F>
    void forEachUsingSameContext(F f) const {
        f(*this);
    }

    static FinegrainedSetupParam makeInactive() {
        FinegrainedSetupParam res{0,0,0,0};
        res.setCost(0.f);
        return res;
    }
};

template<typename T, typename Tag>
struct FinegrainedPartitionnedFFTConvolutionSimulation {
    using SetupParam = FinegrainedSetupParam;
    
    void setup(SetupParam const & p) {
        reset();
        this->p = p;

        partition_size_minus_one = p.partition_size - 1;

        count_multiplicative_grains = p.countMultiplicativeGrains();
        granularity = p.getGranularity();
        
        costXForwardFft = fft::AlgoCosts<Tag,T>::cost_fft_forward(get_fft_length());
        costXInverseFft = fft::AlgoCosts<Tag,T>::cost_fft_inverse(get_fft_length());
        cost_mult = fft::RealFBinsCosts<Tag,T>::cost_multiply(p.partition_size);
        cost_mult_add = fft::RealFBinsCosts<Tag,T>::cost_multiply_add(p.partition_size);
        cost_add_assign = fft::RealSignalCosts<Tag,T>::cost_add_assign(get_fft_length());
    }
    
    auto const & getParam() const {
        return p;
    }
    
    void setMultiplicationGroupLength(int l) {
        auto p2 = p;
        p2.multiplication_group_size = l;
        setup(p2);
    }
    
    int getBlockSize() const {
        return 1+partition_size_minus_one;
    }

    bool isValid() const {
        if(p.multiplication_group_size == 0) {
            return p.partition_count == 0;
        }
        return p.countGrains() <= p.partition_size;
    }
    
    bool isZero() const {
        return p.partition_count == 0;
    }

    void reset() {
        count_multiplicative_grains = 0;
        granularity = 0;
        
        reset_states();
    }
    
    /*
     Dual method of FinegrainedPartitionnedFFTConvolutionSimulation::step()
     */
    double simuStep(XFFtCostFactors const & xFftCostFactors) {
        double cost = 0.;
        if(unlikely(isZero())) {
            return cost;
        }
        ++grain_counter;
        auto g = nextGrain();
        assert(g.first >= 0);
        if(unlikely(g.first == 0)) {
            cost += simuDoGrain(g.second,
                                xFftCostFactors);
            updatePostGrain(g.second);
            grain_counter = 0;
        }
        ++x_progress; // progress may overflow, that's ok because we use only its lower bits.

        return cost;
    }
    
private:
    int32_t grain_counter = 0;
    int32_t grain_number = 0;
    uint32_t x_progress = 0;

    SetupParam p;
    
    int32_t count_multiplicative_grains = 0;
    int32_t granularity = 0;
    uint32_t partition_size_minus_one = -1;

    double costXForwardFft = 0.;
    double costXInverseFft = 0.;
    double cost_mult = 0.;
    double cost_mult_add = 0.;
    double cost_add_assign = 0;
    
    void reset_states() {
        grain_counter = 0;
        grain_number = 0;
        x_progress = 0;
    }
    
    void updatePostGrain(GrainType g) {
        if(g==GrainType::FFT) {
            this->grain_number = 0;
        }
        else {
            ++this->grain_number;
        }
    }
    
    std::pair<int, GrainType> nextGrain() const {
        
        auto const n_mult_grains_remaining = count_multiplicative_grains - grain_number;
        
        if(unlikely(n_mult_grains_remaining < 0)) {
            int distToFFTGrain = partition_size_minus_one - (x_progress & partition_size_minus_one);
            assert(distToFFTGrain >= 0);
            return {distToFFTGrain, GrainType::FFT};
        }
        
        int const dist = granularity - grain_counter;
        
        assert(dist >= 0);
        if(unlikely(n_mult_grains_remaining == 0)) {
            return {dist, GrainType::IFFT};
        }
        else {
            Assert(n_mult_grains_remaining > 0);
            return {dist, GrainType::MultiplicationGroup};
        }
    }
    
    double simuDoGrain(GrainType g,
                       XFFtCostFactors const & xFftCostFactors) const
    {
        double cost {};

        auto const fft_length = get_fft_length();

        if(g == GrainType::FFT) {
            Assert(0 == ((x_progress+1) & partition_size_minus_one)); // make sure 'rythm is good'

            // the fft is done by x_and_ffts
            cost += costXForwardFft * xFftCostFactors.getCostMultiplicator(fft_length);
        }
        else {
            if(g == GrainType::MultiplicationGroup) {
                auto const M = getMultiplicationsGroupMaxSize();
                int offset = M * grain_number;
                assert(offset >= 0);
                assert(offset < p.partition_count);

                int offset_end = std::min(offset + M,
                                          p.partition_count);
                if(offset == 0) {
                    cost += cost_mult;
                    offset = 1;
                }
                int const nRemaining = std::max(0,
                                                offset_end-offset);
                cost += nRemaining * cost_mult_add;
            }
            else {
                Assert(g==GrainType::IFFT);

                cost += costXInverseFft;
                cost += cost_add_assign;
            }
        }
        return cost;
    }
    
    auto get_fft_length() const { return 2 * (1+partition_size_minus_one); }
    auto getMultiplicationsGroupMaxSize() const { return p.multiplication_group_size; }
};
    
template<typename T, typename FFTTag>
struct Simulation_<FinegrainedSetupParam, T, FFTTag> {
    using type = FinegrainedPartitionnedFFTConvolutionSimulation<T, FFTTag>;
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
    static constexpr bool step_can_error = false;
    
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
            int const n_grains = countGrains();
            int const granularity = n_grains ? getBlockSize()/n_grains : 0;

            os << "progress [_/" << getBlockSize()
            << "] grain_counter " << grain_counter << "/" << granularity
            << " grain " << getGrainNumber() << "/" << countGrains() << std::endl;

            doLogComputeState(os);
        }
    }
    
    void setup(SetupParam const & p) {
        set_partition_size(p.partition_size, p.partition_count);
        setMultiplicationGroupLength(p.multiplication_group_size);
    }
    std::array<int, nComputePhaseable> getComputePeriodicities() const {
        return {getBlockSize()};
    }
    // in [0, getComputePeriodicity())
    std::array<int, nComputePhaseable> getComputeProgresses() const {
        auto const sz = getBlockSize();
        auto const res = sz ? static_cast<int>(progress % sz) : 0;
        return {res};
    }
    void setComputeProgresses(std::array<int, nComputePhaseable> const & progresses) {
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
    using Algo = typename Parent::Algo;
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
    using Parent::doFlushToSilence;
    using Parent::get_multiply_add_result;
    
    void setCoefficients(a64::vector<T> coeffs_) {
        reset();
        
        auto const N = coeffs_.size();
        auto const fft_length = get_fft_length(N);
        if(!fft_length) {
            return;
        }
        fft.setContext(Contexts::getInstance().getBySize(fft_length));
        
        result.resize(fft_length);
        x.resize(fft_length);
        y.resize(fft_length);
        
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
        Parent::doFlushToSilence();
        if(!result.empty()) {
            zero_signal(result);
        }
        progress = 0;
        if(!y.empty()) {
            zero_signal(y);
        }
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
        progress = 0;
        if(!y.empty()) {
            zero_signal(y);
        }
        grain_counter = 0;
    }
    
public:
    
    // Just used for calibration
    void fastForwardToComputation(GrainType t, T val = 1) {
        switch(t) {
            case GrainType::FFT:
                while(progress != getBlockSize()-1) {
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
                while(progress != 0) {
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
        return (progress == getBlockSize()-1) || (grain_counter+1 == getBlockSize()/countGrains());
    }
    
    T step(T val) {
        assert(isValid());
        auto g = nextGrain();
        assert(g.first > 0);
        
        x[progress] = typename RealSignal::value_type(val);
        ++progress;
        ++grain_counter;
        
        if(g.first == 1 && g.second) {
            doGrain(*g.second);
        }
        
        assert(progress < getBlockSize());
        assert(progress >= 0);
        assert(!y.empty());
        
        return get_signal(y[progress]);
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
                x[progress] = typename RealSignal::value_type(input_buffer[i]);
                ++progress;
                if(i == g.first-1 && g.second) {
                    doGrain(*g.second);
                }
                f(output_buffer[i], get_signal(y[progress]));
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
                x[progress] = typename RealSignal::value_type(0);
                ++progress;
                if(i == g.first-1 && g.second) {
                    doGrain(*g.second);
                }
                output_buffer[i] += get_signal(y[progress]);
            }
            if(0 == nSamples) {
                return;
            }
            output_buffer += g.first;
        }
    }
    
private:
    
    std::pair<int, std::optional<GrainType>> nextGrain() const {
        int const block_size = getBlockSize();
        int distanceToFFTGrain = block_size - progress;
        assert(distanceToFFTGrain >= 0);
        
        int const n_grains = countGrains();
        auto cur_grain = getGrainNumber();
        if(cur_grain >= n_grains) {
            // spread is not optimal
            return {distanceToFFTGrain, GrainType::FFT};
        }
        
        int const granularity = block_size/n_grains;
        int distanceToOtherGrain = granularity - grain_counter;
        
        // in case of equality, FFT wins.
        if(distanceToOtherGrain < distanceToFFTGrain) {
            assert(distanceToOtherGrain >= 0);
            assert(cur_grain <= n_grains);
            if( cur_grain < n_grains - 1 ) {
                return {distanceToOtherGrain, GrainType::MultiplicationGroup};
            }
            else {
                Assert(cur_grain == n_grains - 1);
                return {distanceToOtherGrain, GrainType::IFFT};
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
                assert(progress == block_size); // make sure 'rythm is good', i.e we exhausted the first half of the y vector
                progress = 0;
                
                //////////////////////// FFT grain /////////////////////////////////////
                
                // x is already padded
                
                compute_x_fft(fft, x);
                
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
                copy(&y[block_size],
                     &result[block_size],
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
                fft.inverse(get_multiply_add_result().data(),
                            result.data(),
                            get_fft_length());
                increment_grain();
                break;
            }
        }
        grain_counter = 0;
    }
    
private:
    int grain_counter = 0;
    int progress = 0;
    Algo fft;
    RealSignal x, y, result;
};
      
constexpr Latency minLatencyLateHandlerWhenEarlyHandlerIsDefaultOptimizedFIRFilter {
    ScaleConvolution_::latencyForDroppedConvolutions(ScaleConvolution_::nDroppedOptimalFor_Split_Bruteforce_Fft)
};
    
template<typename SetupParam>
int constexpr getMinLg2PartitionSz() {
    int partition_sz = 1;
    for(;;partition_sz *= 2) {
        if(SetupParam::getLatencyForPartitionSize(partition_sz) >= minLatencyLateHandlerWhenEarlyHandlerIsDefaultOptimizedFIRFilter) {
            break;
        }
    }
    
    return power_of_two_exponent(partition_sz);
}


template <typename T, template<typename> typename Allocator, typename Tag>
struct FinegrainedPartitionnedFFTConvolutionCRTP {
    using FPT = T;
    using FFTTag = Tag;
    
    using RealSignal = typename fft::RealSignal_<Tag, FPT>::type;
    using Signal_value_type = typename fft::RealSignal_<Tag, FPT>::value_type;
    
    using RealFBins = typename fft::RealFBins_<Tag, FPT, Allocator>;
    using CplxFreqs = typename RealFBins::type;
    static constexpr auto zero = RealFBins::zero;
    
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
    static constexpr Latency getLatencyForPartitionSize(int partSz) {
        return Latency(2*partSz - 1);
    }
    bool handlesCoefficients() const {
        return countPartitions() * getBlockSize() > 0;
    }
    Latency getLatency() const {
        Assert(handlesCoefficients());
        return getLatencyForPartitionSize(partition_size);
    }
    auto getGranularMinPeriod() const { return getBlockSize() / countGrains(); }
    bool isValid() const { return partition_size > 0 && mult_grp_len > 0 && countGrains() <= getBlockSize(); }
    
    int countPartitions() const { return partition_count; }
    
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

    void doFlushToSilence() {
        reset_base_states();
        
        if(!work.empty()) {
            zero(work);
        }

        for(auto & fft_of_delayed_x : ffts_of_delayed_x) {
            zero(fft_of_delayed_x);
        }
        if(!ffts_of_delayed_x.empty()) {
            ffts_of_delayed_x.setByIndex(0);
        }
    }
    
    // 0 means that it will not be used.
    void set_partition_size(int sz, int count) {
        partition_size = sz;
        assert(sz == 0 || is_power_of_two(sz));
        partition_count = count;
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
        auto const fft_length = get_fft_length();

        work.resize(fft_length);

        assert(partition_size > 0);
        
        auto const n_partitions = imajuscule::countPartitions(coeffs_.size(), partition_size);
        if(n_partitions != countPartitions()) {
            throw std::logic_error("inconsistent partition sizes");
        }
        // if one partition is partial, it will be padded with zeroes.
        coeffs_.resize(n_partitions * partition_size);
        
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
    
protected:
    
    int getGrainNumber() const { return grain_number; }
    
    void compute_x_fft(Algo const & fft, RealSignal const & x) {
        assert(grain_number == countGrains());
        grain_number = 0;
        ffts_of_delayed_x.go_back();
        auto & oldest_fft_of_delayed_x = *ffts_of_delayed_x.cycleEnd();
        auto const fft_length = get_fft_length();
        assert(fft_length == oldest_fft_of_delayed_x.size());
        fft.forward(x.begin(), oldest_fft_of_delayed_x.data(), fft_length);
    }
    
    void do_some_multiply_add() {
        auto const M = getMultiplicationsGroupMaxSize();
        auto const offset_base = M * (grain_number - 1);
        assert(offset_base >= 0);
        assert(offset_base < ffts_of_partitionned_h.size());
        auto it_fft_of_partitionned_h = ffts_of_partitionned_h.begin() + offset_base;
        for(auto offset = offset_base, offset_end = std::min((offset_base+M), static_cast<int>(ffts_of_partitionned_h.size()));
            offset != offset_end;
            ++offset, ++it_fft_of_partitionned_h) {
            assert(it_fft_of_partitionned_h < ffts_of_partitionned_h.end());
            
            auto const & fft_of_delayed_x = ffts_of_delayed_x.get_forward(offset);
            
            if(offset == 0) {
                RealFBins::multiply(work.data()                /*   =   */,
                                    fft_of_delayed_x.data(),   /*   x   */   it_fft_of_partitionned_h->data(),
                                    partition_size);
            }
            else {
                RealFBins::multiply_add(work.data()                /*   +=   */,
                                        fft_of_delayed_x.data(),   /*   x   */   it_fft_of_partitionned_h->data(),
                                        partition_size);
            }
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
    int partition_count = 0;
    int grain_number = 0;
    cyclic<CplxFreqs> ffts_of_delayed_x;
    std::vector<CplxFreqs> ffts_of_partitionned_h;
    
    CplxFreqs work;
};

template <typename T, template<typename> typename Allocator, typename FFTTag>
using FinegrainedPartitionnedFFTConvolution =
    FinegrainedFFTConvolutionBase< FinegrainedPartitionnedFFTConvolutionCRTP<T, Allocator, FFTTag> >;

  /*
   input parameters :
   - sz_audio_cb, n_channels
   - impulse response length

   output parameters:
   - lg(partition size)

   1D - Gradient descent according to cost 'max grain computation time' with variable parameters 'lg(partition_size)'
   deducing 'number of multiplications per grain' by finding the parameter that leads to grain computations time just below max(fft, ifft),
   with the constraint that one computation at most occurs per 'equivalent' audio callback.
   */
template<typename Sim, typename SetupParam = typename Sim::SetupParam, typename GradientDescent = GradientDescent<SetupParam>>
auto find_optimal_partition_size(GradientDescent & gradient_descent,
                                 int const n_iterations,
                                 int const n_channels,
                                 int const n_scales,
                                 int const sz_audio_cb,
                                 int const zero_latency_response_size,
                                 std::ostream & os) -> std::optional<SetupParam>
{
    auto n_coeffs_for_latency = [zero_latency_response_size, n_scales](Latency const latency) -> std::optional<int> {
        if( latency < minLatencyLateHandlerWhenEarlyHandlerIsDefaultOptimizedFIRFilter) {
            // invalid case. We pass 'getMinLg2PartitionSz()' to 'run' so that
            // this case doesn't happen on the first try.
            return {};
        }
        // we substract the count of coefficients handled by the early coefficients handler.
        // it favors long partitions because we don't take into account
        // the cost of the early coefficient handler, but for long responses, where optimization matters,
        // the induced bias is negligible.
        
        auto late_response_sz = std::max(0,
                                         zero_latency_response_size - latency.toInteger());
        
        auto scale_sz = SameSizeScales::get_scale_sz(late_response_sz, n_scales);
        if(scale_sz <= 0) {
            return {};
        }
        if(n_scales > 1) {
            // verify that we can have more than one scale (i.e the delay needed to scale is strictly positive)
            // in theory we could relax the constraint (0 delays are ok but the implementation
            // doesn't support that).
            if(SameSizeScales::getDelays(scale_sz, latency) <= 0) {
                return {};
            }
        }
        return {scale_sz};
    };

    if(n_scales != 1) {
        throw std::runtime_error("optimization with scales not implemented");
    }

    gradient_descent.setFunction( [sz_audio_cb, n_coeffs_for_latency, n_iterations, n_channels, &os] (int const lg2_partition_size, auto & val)
    {
        using namespace profiling;
        using namespace std;
        using namespace std::chrono;
        
        if(lg2_partition_size < 0) {
            return ParamState::OutOfRange;
        }
        if(lg2_partition_size > 20) {
            throw logic_error("Gradient descent is not working?");
        }
        int const partition_size = pow2(lg2_partition_size);
        //            cout << "partition size : " << partition_size << endl;
        
        auto maybe_impulse_sz = n_coeffs_for_latency(SetupParam::getLatencyForPartitionSize(partition_size));
        if(!maybe_impulse_sz) {
            return ParamState::OutOfRange;
        }
        int const length_impulse = *maybe_impulse_sz;
        
        // the value for multiplication group size is not very important (it will be overriden later on)
        // but needs to lead to a valid FinegrainedSetupParam. We use the highest valid number:
        int const n_partitions = countPartitions(length_impulse,
                                                 partition_size);
        
        Sim sim;
        sim.setup({
            partition_size,
            n_partitions,
            n_partitions,
            0
        });
        
        Assert(sim.getBlockSize() == partition_size);
        
        // TODO [early coefficients cost] substract the early coefficients from length_impulse
        if(!sim.isValid()) {
            return ParamState::OutOfRange;
        }
        
        cyclic<double> frames_costs{static_cast<size_t>(partition_size)};
        cyclic<double> phased_grains_costs(static_cast<size_t>(partition_size));
        
        struct CostForPhase {
            double getCost() const {
                return cost;
            }
            double cost = 0.;
            double phase = 0.;
        };
        
        RangedGradientDescent<CostForPhase> rgd([sz_audio_cb, n_channels, &frames_costs, &phased_grains_costs, &sim](int multiplication_group_size, auto & result)
        {
            sim.setMultiplicationGroupLength(multiplication_group_size);
            if(!sim.isValid()) {
                return ParamState::OutOfRange;
            }

            XFFtCostFactors emptyCostFactors; // in the future we wil take the number of sources into account here.
            
            for(int i=0, end=frames_costs.size(); i<end; ++i) {
                *(frames_costs.begin() + i) = sim.simuStep(emptyCostFactors);
            }
            result.cost = computeMaxSlidingSum(frames_costs,
                                               sz_audio_cb);
            result.phase = 0.;
            //LG(INFO,"phase=0 half grains cost %f", cost);
            result.cost *= n_channels;
            // now cost is the 'phase == 0' cost

            if(n_channels >= 2) {
                int const nMinFullCbInGranularity = sim.getParam().getGranularity() / sz_audio_cb;
                if(nMinFullCbInGranularity >= n_channels) {
                    // naive phasing leads to an optimal solution : with a phase of 'sz_audio_cb',
                    // during every callback call there will be _at most_ one grain among all channels
                    // that will be computed.
                    result.phase = sz_audio_cb;

                    compute_phased_sum(frames_costs,
                                       result.phase,
                                       n_channels,
                                       phased_grains_costs);
                    result.cost = computeMaxSlidingSum(phased_grains_costs,
                                                       sz_audio_cb);
                }
                else {
                    // naive phasing is not enough to have an optimal solution
                    // because the grains are too close to each other
                    
                    for(int i=0; i<2; ++i) {
                        int base_phase = 0;
                        if(i) {
                            base_phase = sim.getParam().getGranularity() / n_channels;
                            if(base_phase == 0) {
                                // redundant with i==0
                                continue;
                            }
                        }
                        
                        int const maxPhase = 1 + (frames_costs.size() / n_channels);
                        int const phaseIncrement = sim.getParam().getGranularity();
                        for(int phase = base_phase;
                            phase <= maxPhase;
                            phase += phaseIncrement)
                        {
                            if(0 == phase) {
                                continue;
                            }
                            compute_phased_sum(frames_costs,
                                               phase,
                                               n_channels,
                                               phased_grains_costs);
                            
                            auto phased_cost = computeMaxSlidingSum(phased_grains_costs,
                                                                    sz_audio_cb);
                            if(phased_cost < result.cost) {
                                result.cost = phased_cost;
                                result.phase = phase;
                            }
                        }
                    }
                }
            }
            
            // cost is now the sum of costs of each channel
            // but cost should be per sample, not per frame, so
            // we divide by the number of channels
            result.cost /= static_cast<float>(n_channels);

            result.cost /= sz_audio_cb;
            // cost == 'worst computation time over one callback, averaged per sample'

            return ParamState::Ok;
        });
        
        range<int> const multiplication_group_length {
            sim.getParam().getLowestValidMultiplicationsGroupSize(),
            sim.getParam().getHighestValidMultiplicationsGroupSize()
        };
        
        CostForPhase best;
        val.multiplication_group_size = rgd.findLocalMinimum(n_iterations, multiplication_group_length, best);
        val.setCost(best.cost);
        val.setPhase(best.phase);
        val.partition_size = partition_size;
        val.partition_count = countPartitions(length_impulse, partition_size);
        
        constexpr auto debug = false;
        if(debug) {
            rgd.plot(true, os);
            rgd.make_exhaustive(multiplication_group_length, os);
            rgd.plot(true, os);
        }
        return ParamState::Ok;
    });
    
    // must yield a valid result:
    int const min_lg2_partitionsz = getMinLg2PartitionSz<SetupParam>();
    
    Optional<SetupParam> min_val;
    auto res = gradient_descent.findLocalMinimum(n_iterations,
                                                 min_lg2_partitionsz,
                                                 min_val);
    if(res) {
        assert(min_val);
        assert(min_val->partition_size == pow2(*res));
    }
    return min_val;
}
    
template<typename T, typename Tag>
struct PartitionAlgo< FinegrainedSetupParam, T, Tag > {
    using Sim = FinegrainedPartitionnedFFTConvolutionSimulation<T, Tag>;
    using SetupParam = FinegrainedSetupParam;
    
    static std::optional<SetupParam> run(int const n_channels,
                                         int const n_scales,
                                         int const n_audio_frames_per_cb,
                                         int const zero_latency_response_size,
                                         std::ostream & os)
    {
        os << "Optimization of FinegrainedPartitionnedFFTConvolution for " << n_scales << " scale(s)" << std::endl;
        IndentingOStreambuf i(os);
        
        assert(n_channels > 0);
        GradientDescent<SetupParam> gd;
        constexpr auto n_iterations = 1;
        std::optional<SetupParam> res = find_optimal_partition_size<Sim>(gd,
                                                                         n_iterations,
                                                                         n_channels,
                                                                         n_scales,
                                                                         n_audio_frames_per_cb,
                                                                         zero_latency_response_size,
                                                                         os);
        constexpr auto debug_gradient_descent = false;
        if constexpr (debug_gradient_descent) {
            os << "Gradient descent report :" << std::endl;
            gd.debug(true, os); // this takes time because it makes the search exhaustive
        }
        return res;
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
