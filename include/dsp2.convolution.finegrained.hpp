
namespace imajuscule {


struct FinegrainedSetupParam2 : public Cost {
    FinegrainedSetupParam2(int partitionSz,
                           int partition_count,
                           int multiplication_group_size,
                           int phase)
    : Cost(phase)
    , multiplication_group_size(multiplication_group_size)
    , partition_size(partitionSz)
    , partition_count(partition_count)
    {}
    
    int multiplication_group_size = 0;
    int partition_size = 0;
    int partition_count = 0;
    
    void logSubReport(std::ostream & os) const override {
        os << "Finegrained, " << partition_count << " partitions of size " << partition_size << ", mult group: " << multiplication_group_size << std::endl;
    }
};

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
    using Tag = typename Parent::FFTTag;
    using Desc = DescFinegrainedFFTConvolutionBase;
    using Algo = AlgoFinegrainedFFTConvolutionBase<typename Parent::Algo>;
    
    using RealSignal = typename fft::RealSignal_<Tag, FPT>::type;
    
    static constexpr auto zero_signal = fft::RealSignal_<Tag, FPT>::zero;

    using Parent::clear;
    using Parent::doLogComputeState;
    using Parent::doSetCoefficients;

    void logComputeState(Algo const & algo, std::ostream & os) const {
        os << "Finegrained ";
        if(isZero()) {
            os << "zero" << std::endl;
        }
        else {
            os << "grain {" << grain_counter << "/" << algo.countGrains()
            << "} progress [" << this->progress << "/" << algo.getBlockSize() << "]" << std::endl;
            doLogComputeState(algo, os);
        }
    }
    
    MinSizeRequirement setCoefficients(Algo const & algo, a64::vector<T> coeffs_) {        
        reset();
        result.resize(algo.get_fft_length());
        
        return doSetCoefficients(algo, std::move(coeffs_));
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
    }
    
    bool isZero() const {
        return result.empty();
    }
    
    int grain_counter = 0;
    int progress = 0;
    RealSignal result;
    
private:
    void reset_states() {
        Parent::reset_base_states();
        progress = 0;
        grain_counter = 0;
    }
};

template <typename Parent>
struct AlgoFinegrainedFFTConvolutionBase : public Parent {
    using State = StateFinegrainedFFTConvolutionBase<typename Parent::State>;
    
    using T = typename Parent::FPT;
    using FPT = T;
    using Tag = typename Parent::FFTTag;
    using Desc = DescFinegrainedFFTConvolutionBase;
    
    using SetupParam = FinegrainedSetupParam2;
    
    static constexpr auto add_assign = fft::RealSignal_<Tag, FPT>::add_assign;
    
    using RealSignal = typename fft::RealSignal_<Tag, FPT>::type;
    
    using Parent::countPartitions;
    using Parent::countGrains;
    using Parent::do_some_multiply_add;
    using Parent::get_multiply_add_result;
    using Parent::get_fft_length;
    using Parent::getBlockSize;
    using Parent::getGranularMinPeriod;
    using Parent::getLatencyForPartitionSize;
    using Parent::isValid;
    using Parent::set_partition_info;
    using Parent::setMultiplicationGroupLength;
    
    void setup(SetupParam const & p) {
        set_partition_info(p.partition_size, p.partition_count);
        setMultiplicationGroupLength(p.multiplication_group_size);
    }
    
    void step(State & s,
              XAndFFTS<T, Tag> const & x_and_ffts,
              Y<T, Tag> & y) {
        if(unlikely(countPartitions() == 0)) {
            return;
        }
        assert(isValid());
        auto g = nextGrain(s);
        assert(g.first > 0);
        
        ++s.progress;
        ++s.grain_counter;
        
        if(g.first == 1 && g.second) {
            doGrain(s, x_and_ffts, y, *g.second);
        }
        
        assert(s.progress < getBlockSize());
        assert(s.progress >= 0);
    }
    
private:
    
    std::pair<int, std::optional<GrainType>> nextGrain(State const & s) const {
        int const block_size = getBlockSize();
        int distanceToFFTGrain = block_size - s.progress;
        assert(distanceToFFTGrain >= 0);
        
        int const n_grains = countGrains();
        auto cur_grain = s.grain_number;
        if(cur_grain >= n_grains) {
            return {distanceToFFTGrain, GrainType::FFT};
        }
        
        int const granularity = block_size/n_grains;
        int distanceToOtherGrain = granularity - s.grain_counter;
        
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
    
    void doGrain(State & s,
                 XAndFFTS<T, Tag> const & x_and_ffts,
                 Y<T, Tag> & y,
                 GrainType g)
    {
        switch(g)
        {
            case GrainType::FFT:
            {
                // the fft is implicit (in x_and_ffts)
                
                int const N = getBlockSize();
                assert(s.progress == N); // make sure 'rythm is good', i.e we exhausted the first half of the y vector
                
                Assert(y.uProgress+(N-1) < y.y.size());
                add_assign(y.y.begin() + y.uProgress,
                           s.result.begin(),
                           N);
                
                s.progress = 0;
                s.grain_number = 1;
                break;
            }
            case GrainType::MultiplicationGroup:
            {
                auto const fft_length = get_fft_length();
                auto const & ffts = x_and_ffts.find_ffts(fft_length);
                
                do_some_multiply_add(s, ffts);
                ++s.grain_number;
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
                
                ffts.fft.inverse(get_multiply_add_result(), s.result, fft_length);
                
                ++s.grain_number;
                break;
            }
            case GrainType::Nothing:
                Assert(0);
                break;
        }
        s.grain_counter = 0;
    }
};


/*
 */


template <typename T, typename Tag>
struct AlgoFinegrainedPartitionnedFFTConvolutionCRTP;

template <typename T, typename Tag>
struct StateFinegrainedPartitionnedFFTConvolutionCRTP {
    using Algo = AlgoFinegrainedPartitionnedFFTConvolutionCRTP<T, Tag>;
    
    using FPT = T;
    using FFTTag = Tag;
    
    using RealSignal = typename fft::RealSignal_<Tag, FPT>::type;
    using Signal_value_type = typename fft::RealSignal_<Tag, FPT>::value_type;
    
    static constexpr auto scale = fft::RealFBins_<Tag, FPT>::scale;
    static constexpr auto multiply = fft::RealFBins_<Tag, FPT>::multiply;
    static constexpr auto multiply_add = fft::RealFBins_<Tag, FPT>::multiply_add;
    
    using CplxFreqs = typename fft::RealFBins_<Tag, FPT>::type;
    
    void doLogComputeState(Algo const & algo, std::ostream & os) const {
        os << algo.countPartitions() << " partitions "
        << algo.getMultiplicationsGroupMaxSize() << " multGroupMaxSize" << std::endl;
    }
    
    bool empty() const {
        return ffts_of_partitionned_h.empty();
    }
    void clear() {
        ffts_of_partitionned_h.clear();
    }
    
    double getEpsilon(Algo const & algo) const {
        return algo.countPartitions() * (fft::getFFTEpsilon<FPT>(algo.get_fft_length()) + 2 * std::numeric_limits<FPT>::epsilon());
    }
    
public:
    MinSizeRequirement doSetCoefficients(Algo const & algo,
                                         a64::vector<T> coeffs_) {
        auto const N = algo.getBlockSize();
        assert(N > 0);
        
        auto const n_partitions = imajuscule::countPartitions(coeffs_.size(), N);
        if(n_partitions != algo.countPartitions()) {
            throw std::logic_error("wrong number of partitions");
        }
        // if one partition is partial, it will be padded with zeroes.
        coeffs_.resize(n_partitions * N);
        
        ffts_of_partitionned_h.resize(n_partitions);
        
        auto const fft_length = algo.get_fft_length();
        
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
                fft.forward(coeffs_slice.begin(), fft_of_partitionned_h, fft_length);
                scale(fft_of_partitionned_h, factor);
            }
        }
        Assert(it_coeffs == coeffs_.end());
        
        Assert(fft_length > 1);
        return {
            static_cast<int>(fft_length/2), // x block size
            static_cast<int>(fft_length/2), // y block size
            static_cast<int>(fft_length/2), // y anticipated writes (because we write "in the future" of y in the ifft step)
            {
                {fft_length, n_partitions}
            }
        };
    }

protected:
    void reset_base_states() {
        grain_number = 1;
    }
    
public:
    int grain_number = 0;
    std::vector<CplxFreqs> ffts_of_partitionned_h;
};


template <typename T, typename Tag>
struct AlgoFinegrainedPartitionnedFFTConvolutionCRTP {
    using State = StateFinegrainedPartitionnedFFTConvolutionCRTP<T, Tag>;
    
    using FPT = T;
    using FFTTag = Tag;
    
    using RealSignal = typename fft::RealSignal_<Tag, FPT>::type;
    using Signal_value_type = typename fft::RealSignal_<Tag, FPT>::value_type;
    
    using CplxFreqs = typename fft::RealFBins_<Tag, FPT>::type;

    static constexpr auto multiply = fft::RealFBins_<Tag, FPT>::multiply;
    static constexpr auto multiply_add = fft::RealFBins_<Tag, FPT>::multiply_add;
    
    auto get_fft_length() const { return 2 * partition_size; }
    
    int getWriteYBlockSize() const { return partition_size; }
    auto getBlockSize() const { return partition_size; }
    static constexpr int getLatencyForPartitionSize(int partSz) {
        return 2*partSz - 1;
    }
    auto getLatency() const { return getLatencyForPartitionSize(partition_size); }
    auto getGranularMinPeriod() const { return getBlockSize() / countGrains(); }
    bool isValid() const { return partition_size > 0 && mult_grp_len > 0 && countGrains() <= getBlockSize(); }
     
    int countPartitions() const { return partition_count; }
    
    double getEpsilon() const {
        return countPartitions() * (fft::getFFTEpsilon<FPT>(get_fft_length()) + 2 * std::numeric_limits<FPT>::epsilon());
    }
    
protected:
    void setMultiplicationGroupLength(int length) {
        mult_grp_len = length;
    }
    
    // sz = 0 or count = 0 means that it will not be used.
    void set_partition_info(int sz, int count) {
        partition_count = count;
        partition_size = sz;
        assert(sz == 0 || is_power_of_two(sz));
        
        work.resize(get_fft_length());
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
            if( (countPartitions() - 1)/i <= diff) {
                return i;
            }
        }
    }
    
    int getHighestValidMultiplicationsGroupSize() const { return countPartitions(); }
    
protected:
    void do_some_multiply_add(State & s,
                              typename XAndFFTS<T, Tag>::FFTs const & ffts) {
        auto const M = getMultiplicationsGroupMaxSize();
        auto const offset_base = M * (s.grain_number - 1);
        assert(offset_base >= 0);
        assert(offset_base < s.ffts_of_partitionned_h.size());
        auto it_fft_of_partitionned_h = s.ffts_of_partitionned_h.begin() + offset_base;
        for(auto offset = offset_base, offset_end = std::min((offset_base+M), static_cast<int>(s.ffts_of_partitionned_h.size()));
            offset != offset_end;
            ++offset, ++it_fft_of_partitionned_h)
        {
            assert(it_fft_of_partitionned_h < s.ffts_of_partitionned_h.end());
            
            auto const & fft_of_delayed_x = ffts.ffts.get_backward(offset);
            
            if(offset == 0) {
                multiply(work                /*   =   */,
                         fft_of_delayed_x,   /*   x   */   *it_fft_of_partitionned_h);
            }
            else {
                multiply_add(work                /*   +=   */,
                             fft_of_delayed_x,   /*   x   */   *it_fft_of_partitionned_h);
            }
        }
    }
    
    auto const & get_multiply_add_result() const {
        return work;
    }
    
private:
    int mult_grp_len = 0;
    int partition_size = -1;
    int partition_count = 0;
    
    CplxFreqs work;
};

template <typename T, typename FFTTag>
using AlgoFinegrainedPartitionnedFFTConvolution = AlgoFinegrainedFFTConvolutionBase< AlgoFinegrainedPartitionnedFFTConvolutionCRTP<T, FFTTag> >;

}
