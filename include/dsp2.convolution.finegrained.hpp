
namespace imajuscule {

struct DescFinegrainedFFTConvolutionBase {
    static constexpr int nCoefficientsFadeIn = 0;
    static constexpr bool has_subsampling = false;
    static constexpr bool step_can_error = false;
};

template <typename Parent>
struct AlgoFinegrainedFFTConvolutionBase;

template <typename Parent>
struct StateFinegrainedFFTConvolutionBase : public Parent {
    using T = typename Parent::FPT;
    using FPT = T;
    using Tag = typename Parent::Tag;
    using Desc = DescFinegrainedFFTConvolutionBase;
    using Algo = AlgoFinegrainedFFTConvolutionBase<typename Parent::Algo>;
    
    using RealSignal = typename fft::RealSignal_<Tag, FPT>::type;
    
    static constexpr auto zero_signal = fft::RealSignal_<Tag, FPT>::zero;

    using Parent::clear;
    using Parent::doSetCoefficients;
    using Parent::doLogComputeState;
    using Parent::doOnContextFronteer;
    using Parent::doFlushToSilence;
    
    void logComputeState(Algo const & algo, std::ostream & os) const {
        os << "Finegrained ";
        if(isZero()) {
            os << "zero" << std::endl;
        }
        else {
            os << "progress [_/" << algo.getBlockSize()
            << "] grain_counter " << grain_counter << "/" << algo.getGranularity();
            doLogComputeState(algo, os);
        }
    }
    
    void setCoefficients(Algo const & algo, a64::vector<T> coeffs_) {
        reset();
        result.resize(algo.get_fft_length());
        
        doSetCoefficients(algo, std::move(coeffs_));
    }
    
    template<typename F>
    void onContextFronteer(F f) {
        doOnContextFronteer(f);
    }
    
    void reset() {
        reset_states();
        result.clear();
        clear();
    }
    void flushToSilence() {
        reset_states();

        if(!result.empty()) {
            zero_signal(result);
        }
        
        doFlushToSilence();
    }
    
    bool isZero() const {
        return result.empty();
    }
    
    int grain_counter = 0;
    mutable RealSignal result;
    
private:
    void reset_states() {
        Parent::reset_base_states();
        grain_counter = 0;
    }
};

template <typename Parent>
struct AlgoFinegrainedFFTConvolutionBase : public Parent {
    using State = StateFinegrainedFFTConvolutionBase<typename Parent::State>;
    
    using T = typename Parent::FPT;
    using FPT = T;
    using Tag = typename Parent::Tag;
    using Desc = DescFinegrainedFFTConvolutionBase;
    
    static constexpr auto add_assign = fft::RealSignal_<Tag, FPT>::add_assign;
    
    using RealSignal = typename fft::RealSignal_<Tag, FPT>::type;
    
    template<typename TT>
    using Allocator = typename Parent::template Allocator<TT>;
    
    using Parent::countPartitions;
    using Parent::countGrains;
    using Parent::do_some_multiply_add;
    using Parent::get_fft_length;
    using Parent::getBlockSize;
    using Parent::isValid;
    using Parent::setMultiplicationGroupLength;
    using Parent::getGranularity;
    
    void dephaseSteps(State & s,
                      int const n_steps) const {
        //LG(INFO, "dephase by %d", n_steps);
        if(s.isZero()) {
            return;
        }
        if(unlikely(countPartitions() == 0)) {
            return;
        }
        for(int x_progress = 0;
            x_progress < n_steps;
            ++x_progress)
        {
            // Note that the high bits of x_progress will be ignored (& mask with block_size-1)
            auto g = nextGrain(s, x_progress);
            if(g.first == 1 && g.second) {
                s.updatePostGrain(*g.second);
                s.grain_counter = 0;
            }
            else {
                ++s.grain_counter;
            }
        }
    }

    template<template<typename> typename Allocator2>
    void step(State & s,
              XAndFFTS<T, Allocator2, Tag> const & x_and_ffts,
              Y<T, Tag> & y) const {
        if(s.isZero()) {
            return;
        }
        if(unlikely(countPartitions() == 0)) {
            return;
        }
        assert(isValid());
        auto g = nextGrain(s, x_and_ffts.progress);
        assert(g.first > 0);

        if(g.first == 1 && g.second) {
            doGrain(s, x_and_ffts, y, *g.second);
            s.updatePostGrain(*g.second);
            s.grain_counter = 0;
        }
        else {
            ++s.grain_counter;
        }
    }

    void flushToSilence(State & s) const {
        s.flushToSilence();
    }

private:
    
    std::pair<int, std::optional<GrainType>> nextGrain(State const & s,
                                                       int const x_progress) const {
        unsigned int const block_size = getBlockSize();
        int distanceToFFTGrain = block_size - ((x_progress-1) & (block_size-1));
        assert(distanceToFFTGrain >= 0);
        
        int const n_grains = countGrains();
        auto cur_grain = s.grain_number;
        if(unlikely(cur_grain >= n_grains)) {
            // spread is not optimal
            return {distanceToFFTGrain, GrainType::FFT};
        }
        
        int distanceToOtherGrain = getGranularity() - s.grain_counter;
        
        // in case of equality, FFT wins.
        if(likely(distanceToOtherGrain < distanceToFFTGrain)) {
            assert(distanceToOtherGrain >= 0);
            assert(cur_grain <= n_grains);
            if(likely(cur_grain < n_grains - 1)) {
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
    
    template<template<typename> typename Allocator2>
    void doGrain(State & s,
                 XAndFFTS<T, Allocator2, Tag> const & x_and_ffts,
                 Y<T, Tag> & y,
                 GrainType g) const
    {
        switch(g)
        {
            case GrainType::FFT:
            {
                // the fft is implicit (in x_and_ffts)
                
                int const N = getBlockSize();
                Assert(0 == (x_and_ffts.progress & (N-1))); // make sure 'rythm is good', i.e we exhausted the first half of the y vector
                
                Assert(y.uProgress+(N-1) < y.y.size());
                add_assign(y.y.begin() + y.uProgress,
                           s.result.begin(),
                           N);
                break;
            }
            case GrainType::MultiplicationGroup:
            {
                auto const fft_length = get_fft_length();
                auto const & ffts = x_and_ffts.find_ffts(fft_length);
                
                do_some_multiply_add(s, ffts);
                break;
            }
            case GrainType::IFFT:
            {
                {
                    unsigned int const N = getBlockSize();
                    
                    Assert(is_power_of_two(N));
                    
                    unsigned int yFutureLocation = (y.uProgress + N - 1) & ~(N-1);
                    Assert(yFutureLocation <= y.y.size());
                    if(yFutureLocation == y.y.size()) {
                        yFutureLocation = 0;
                    }
                    Assert(yFutureLocation + N <= y.zeroed_up_to);
                    add_assign(y.y.begin() + yFutureLocation,
                               s.result.begin() + N,
                               N);
                }
                
                auto const fft_length = get_fft_length();
                
                auto const & ffts = x_and_ffts.find_ffts(fft_length);
                
                ffts.fft.inverse(s.multiply_add_result.data(), s.result, fft_length);
                break;
            }
        }
    }
};

/*
 */


template <typename T, template<typename> typename Allocator, typename FFTTag>
struct AlgoFinegrainedPartitionnedFFTConvolutionCRTP;

template <typename T, template<typename> typename Allocator, typename FFTTag>
struct StateFinegrainedPartitionnedFFTConvolutionCRTP {
    using Algo = AlgoFinegrainedPartitionnedFFTConvolutionCRTP<T, Allocator, FFTTag>;
    
    using FPT = T;
    using Tag = FFTTag;
    
    using RealSignal = typename fft::RealSignal_<Tag, FPT>::type;
    using Signal_value_type = typename fft::RealSignal_<Tag, FPT>::value_type;
    
    static constexpr auto scale = fft::RealFBins_<Tag, FPT, Allocator>::scale;
    static constexpr auto zero = fft::RealFBins_<Tag, FPT, Allocator>::zero;

    using CplxFreqs = typename fft::RealFBins_<Tag, FPT, Allocator>::type;
    
    bool empty() const {
        return ffts_of_partitionned_h.empty();
    }
    void clear() {
        ffts_of_partitionned_h.clear();
    }
    
    double getEpsilon(Algo const & algo) const {
        return algo.countPartitions() * (fft::getFFTEpsilon<FPT>(algo.get_fft_length()) + 2 * std::numeric_limits<FPT>::epsilon());
    }
    void doLogComputeState(Algo const & algo, std::ostream & os) const {
        os << " grain " << grain_number << "/" << algo.countGrains() << std::endl;
        os << algo.countPartitions() << " partitions "
        << algo.getMultiplicationsGroupMaxSize() << " multGroupMaxSize" << std::endl;
    }

public:
    void doSetCoefficients(Algo const & algo,
                           a64::vector<T> coeffs_)
    {
        auto const fft_length = algo.get_fft_length();
        multiply_add_result.resize(fft_length);

        if(fft_length) {
            auto const N = algo.getBlockSize();
            assert(N > 0);
            
            auto const n_partitions = imajuscule::countPartitions(coeffs_.size(), N);
            if(n_partitions != algo.countPartitions()) {
                throw std::logic_error("wrong number of partitions");
            }
            // if one partition is partial, it will be padded with zeroes.
            coeffs_.resize(n_partitions * N);
            
            ffts_of_partitionned_h.resize(n_partitions);
            
            for(auto & fft_of_partitionned_h : ffts_of_partitionned_h) {
                fft_of_partitionned_h.resize(fft_length);
            }
            
            // compute fft of padded impulse response
            
            auto it_coeffs = coeffs_.begin();
            {
                using FFTAlgo = typename fft::Algo_<Tag, FPT>;
                using Contexts = fft::Contexts_<Tag, FPT>;
                FFTAlgo fft(Contexts::getInstance().getBySize(fft_length));
                
                auto const factor = scaleFactor<FFTAlgo>(static_cast<FPT>(fft_length));
                RealSignal coeffs_slice(fft_length, Signal_value_type(0)); // initialize with zeros (second half is padding)
                for(auto & fft_of_partitionned_h : ffts_of_partitionned_h) {
                    auto end_coeffs = it_coeffs + N;
                    assert(end_coeffs <= coeffs_.end());
                    auto slice_it = coeffs_slice.begin();
                    for(;it_coeffs != end_coeffs; ++it_coeffs, ++slice_it) {
                        using RealT = typename RealSignal::value_type;
                        *slice_it = RealT(*it_coeffs);
                    }
                    
                    // coeffs_slice is padded with 0, because it is bigger than partition_size
                    // and initialized with zeros.
                    fft.forward(coeffs_slice.begin(), fft_of_partitionned_h.data(), fft_length);
                    scale(fft_of_partitionned_h, factor);
                }
            }
            Assert(it_coeffs == coeffs_.end());
            
            Assert(fft_length > 1);
        }
        else {
            ffts_of_partitionned_h.clear();
        }
    }
    
    void doFlushToSilence() {
        if(!multiply_add_result.empty()) {
            zero(multiply_add_result);
        }
    }

    template<typename F>
    void doOnContextFronteer(F f) {
    }
    
    void updatePostGrain(GrainType g) {
        if(g==GrainType::FFT) {
            this->grain_number = 1;
        }
        else {
            ++this->grain_number;
        }
    }
    
protected:
    void reset_base_states() {
        grain_number = 1;
    }
    
public:
    int grain_number = 0;
    std::vector<CplxFreqs> ffts_of_partitionned_h;
    CplxFreqs multiply_add_result; // needs to be in the state, not in algo because the results will accumulate over several steps
};


template <typename T, template<typename> typename Alloc, typename FFTTag>
struct AlgoFinegrainedPartitionnedFFTConvolutionCRTP {
    template<typename TT>
    using Allocator = Alloc<TT>;
    
    using State = StateFinegrainedPartitionnedFFTConvolutionCRTP<T, Alloc, FFTTag>;
    
    using FPT = T;
    using Tag = FFTTag;

    using SetupParam = FinegrainedSetupParam;

    using RealSignal = typename fft::RealSignal_<Tag, FPT>::type;
    using Signal_value_type = typename fft::RealSignal_<Tag, FPT>::value_type;
    
    using CplxFreqs = typename fft::RealFBins_<Tag, FPT, Allocator>::type;

    void setup(SetupParam const & p) {
        assert(p.partition_size == 0 || is_power_of_two(p.partition_size));

        // sz = 0 or count = 0 means that it will not be used.
        partition_count = p.partition_count;
        partition_size = p.partition_size;

        setMultiplicationGroupLength(p.multiplication_group_size);
    }
    auto get_fft_length() const { return 2 * partition_size; }
    
    auto getBlockSize() const { return partition_size; }
    Latency getLatency() const {
        Assert(handlesCoefficients());
        return FinegrainedSetupParam::getLatencyForPartitionSize(partition_size);
    }
    bool isValid() const {
        if(mult_grp_len == 0) {
            return partition_count == 0;
        }
        return countGrains() <= getBlockSize();
    }
    bool handlesCoefficients() const {
        return partition_count > 0;
    }

    int countPartitions() const { return partition_count; }
    
    double getEpsilon() const {
        return countPartitions() * (fft::getFFTEpsilon<FPT>(get_fft_length()) + 2 * std::numeric_limits<FPT>::epsilon());
    }
    
protected:
    void setMultiplicationGroupLength(int length) {
        mult_grp_len = length;
        updateGranularity();
    }
        
public:
    auto getMultiplicationsGroupMaxSize() const { return mult_grp_len; }
    auto countMultiplicativeGrains() const { return 1 + (countPartitions()-1)/getMultiplicationsGroupMaxSize(); }
    static constexpr auto countNonMultiplicativeGrains() { return 2; }
    auto countGrains() const {
        return countNonMultiplicativeGrains() + countMultiplicativeGrains();
    }
    int getGranularity() const { return granularity; }
    
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
        if(countPartitions() == 0) {
            return 0;
        }
        for(int i=1;; ++i) {
            if( (countPartitions() - 1)/i <= diff) {
                return i;
            }
        }
    }
    
    int getHighestValidMultiplicationsGroupSize() const { return countPartitions(); }
    
protected:
    template<template<typename> typename Allocator2>
    void do_some_multiply_add(State & s,
                              FFTs<T, Allocator2, Tag> const & ffts) const {
        auto const M = getMultiplicationsGroupMaxSize();
        auto const offset_base = M * (s.grain_number - 1);
        assert(offset_base >= 0);
        assert(offset_base < s.ffts_of_partitionned_h.size());

        int offset_end = std::min((offset_base+M), static_cast<int>(s.ffts_of_partitionned_h.size()));
        int offset = offset_base;
        if(offset == 0) {
            fft::RealFBins_<Tag, FPT, Allocator>::multiply(s.multiply_add_result.data() /* = */,
                                                           ffts.get_by_age(0), /* x */ s.ffts_of_partitionned_h[0].data(),
                                                           partition_size);
            offset = 1;
        }
        for(; offset != offset_end; ++offset)
        {
            fft::RealFBins_<Tag, FPT, Allocator>::multiply_add(s.multiply_add_result.data() /* += */,
                                                               ffts.get_by_age(offset), /* x  */ s.ffts_of_partitionned_h[offset].data(),
                                                               partition_size);
        }
    }
    
    
private:
    int mult_grp_len = 0;
    int partition_size = -1;
    int partition_count = 0;
    int granularity = 0;
    
    void updateGranularity()
    {
        if(0==mult_grp_len) {
            granularity = 0;
        }
        else if(int n_grains = countGrains()) {
            granularity = getBlockSize()/n_grains;
        }
        else {
            granularity = 0;
        }
    }
};

template <typename T, template<typename> typename Allocator, typename FFTTag>
using AlgoFinegrainedPartitionnedFFTConvolution =
AlgoFinegrainedFFTConvolutionBase< AlgoFinegrainedPartitionnedFFTConvolutionCRTP<T, Allocator, FFTTag> >;

}
