
namespace imajuscule {

struct DescFinegrainedFFTConvolution {
    static constexpr int nCoefficientsFadeIn = 0;
    static constexpr bool has_subsampling = false;
    static constexpr bool step_can_error = false;
};

template <typename T, template<typename> typename Allocator, typename FFTTag>
struct AlgoFinegrainedPartitionnedFFTConvolution;

template <typename T, template<typename> typename Allocator, typename FFTTag>
struct StateFinegrainedPartitionnedFFTConvolution {
    using FPT = T;
    using Tag = FFTTag;
    using Desc = DescFinegrainedFFTConvolution;
    using Algo = AlgoFinegrainedPartitionnedFFTConvolution<T, Allocator, FFTTag>;
    
    using RealSignal = typename fft::RealSignal_<Tag, FPT>::type;
    using Signal_value_type = typename fft::RealSignal_<Tag, FPT>::value_type;
    
    using CplxFreqs = typename fft::RealFBins_<Tag, FPT, Allocator>::type;

    static constexpr auto scale = fft::RealFBins_<Tag, FPT, Allocator>::scale;
    static constexpr auto zero = fft::RealFBins_<Tag, FPT, Allocator>::zero;

    static constexpr auto zero_signal = fft::RealSignal_<Tag, FPT>::zero;

    static int getAllocationSz_SetCoefficients(typename Algo::SetupParam const & p) {
        return p.partition_count ? (p.get_fft_length() * (1 + p.partition_count)) : 0;
    }

    void logComputeState(Algo const & algo, std::ostream & os) const {
        os << "Finegrained ";
        if(isZero()) {
            os << "zero" << std::endl;
        }
        else {
            os << "progress [_/" << algo.getBlockSize()
            << "] grain_counter " << grain_counter << "/" << algo.getGranularity();

            os << " grain " << grain_number << "/" << algo.countGrains() << std::endl;
            os << algo.countPartitions() << " partitions "
            << algo.getMultiplicationsGroupMaxSize() << " multGroupMaxSize" << std::endl;
        }
    }
    
    void setCoefficients(Algo const & algo, a64::vector<T> coeffs_) {
        reset();
        
        ffts_of_partitionned_h.clear();
        multiply_add_result.clear();
        
        auto const fft_length = algo.get_fft_length();
        if(!fft_length) {
            return;
        }
        multiply_add_result.resize(fft_length);
        auto const N = algo.getBlockSize();
        assert(N > 0);
        
        auto const n_partitions = imajuscule::countPartitions(coeffs_.size(), N);
        if(n_partitions != algo.countPartitions()) {
            throw std::logic_error("wrong number of partitions");
        }
        // if one partition is partial, pad with zeroes.
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
    
    template<typename F>
    void onContextFronteer(F f) {
    }
    
    void reset() {
        reset_states();
        ffts_of_partitionned_h.clear();
    }
    
    void flushToSilence() {
        reset_states();

        if(!multiply_add_result.empty()) {
            zero(multiply_add_result);
        }
    }
    
    void updatePostGrain(GrainType g) {
        if(g==GrainType::FFT) {
            this->grain_number = 0;
        }
        else {
            ++this->grain_number;
        }
    }
    
    bool isZero() const {
        return ffts_of_partitionned_h.empty();
    }
    
    double getEpsilon(Algo const & algo) const {
        return algo.countPartitions() * (fft::getFFTEpsilon<FPT>(algo.get_fft_length()) + 2 * std::numeric_limits<FPT>::epsilon());
    }
    
    int32_t grain_counter = 0;
    int32_t grain_number = 0;
    std::vector<CplxFreqs> ffts_of_partitionned_h;
    CplxFreqs multiply_add_result; // the results accumulate over several steps
    
private:
    void reset_states() {
        grain_number = 0;
        grain_counter = 0;
    }
};

template <typename T, template<typename> typename Alloc, typename FFTTag>
struct AlgoFinegrainedPartitionnedFFTConvolution {
    using State = StateFinegrainedPartitionnedFFTConvolution<T, Alloc, FFTTag>;
    
    using FPT = T;
    using Tag = FFTTag;
    using Desc = DescFinegrainedFFTConvolution;
        
    using SetupParam = FinegrainedSetupParam;
    
    template<typename TT>
    using Allocator = Alloc<TT>;
    
    using RealSignal = typename fft::RealSignal_<Tag, FPT>::type;
    
    void setup(SetupParam const & p) {
        assert(p.partition_size == 0 || is_power_of_two(p.partition_size));

        // sz = 0 or count = 0 means that it will not be used.
        partition_count = p.partition_count;
        partition_size_minus_one = p.partition_size - 1;

        mult_grp_len = p.multiplication_group_size;
        count_multiplicative_grains = p.countMultiplicativeGrains();
        granularity = p.getGranularity();
    }
    auto get_fft_length() const { return 2 * (1+partition_size_minus_one); }
    
    auto getBlockSize() const { return 1+partition_size_minus_one; }
    int getBiggestScale() const {
        return getBlockSize();
    }
    
    auto getMultiplicationsGroupMaxSize() const { return mult_grp_len; }
    auto countMultiplicativeGrains() const {
        return count_multiplicative_grains;
    }
    
    auto countGrains() const {
        return count_multiplicative_grains + SetupParam::countNonMultiplicativeGrains();
    }
    int getGranularity() const { return granularity; }
    
    Latency getLatency() const {
        Assert(handlesCoefficients());
        return FinegrainedSetupParam::getLatencyForPartitionSize(getBlockSize());
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
    
    void dephaseStep(State & s,
                     unsigned int const x_progress) const {
        //LG(INFO, "dephase by %d", n_steps);
        if(unlikely(countPartitions() == 0)) {
            return;
        }
        ++s.grain_counter;
        // Note that the high bits of x_progress will be ignored (& mask with block_size-1)
        auto g = nextGrain(s, x_progress);
        if(g.first == 0) {
            s.updatePostGrain(g.second);
            s.grain_counter = 0;
        }
    }

    
    template<template<typename> typename Allocator2, typename WorkData>
    void stepVectorized(State & s,
                        XAndFFTS<T, Allocator2, Tag> const & x_and_ffts,
                        Y<T, Tag> & y,
                        WorkData * workData,
                        int vectorSz) const
    {
        if(unlikely(countPartitions() == 0)) {
            return;
        }
        assert(isValid());

        int const N = getBlockSize();
        
        int const log2Blocksize = count_trailing_zeroes(N);
        int const h = static_cast<unsigned int>(x_and_ffts.progress+1) & ~(N-1);
        
        int yPhase = 0;
        //int const yPhase = (N - ((static_cast<unsigned int>(x_and_ffts.progress+2-vectorSz)) & (N-1))) & (N-1);

        unsigned int x_progress_simu = x_and_ffts.progress + 1 - vectorSz;

        int const l = static_cast<unsigned int>(x_progress_simu) & ~(N-1);
        int nForceSteps = (h-l) >> log2Blocksize;

        do {
            Assert(vectorSz > 0);
                        
            /*
            int nForceSteps2 = x_and_ffts.countFftsSinceForhalfSize(vectorSz,
                                                                   N);
            Assert(nForceSteps==nForceSteps2);
             */
            
            ++s.grain_counter;
            auto g = nextGrain(s, x_progress_simu);
            x_progress_simu += g.first + 1;
            assert(g.first >= 0);
            if(g.first <= vectorSz-1) {
                yPhase += g.first;
                vectorSz -= g.first + 1;
#ifndef NDEBUG
                s.grain_counter += g.first; // just to satisfy an assert, but since we don't use grain_counter there is no need to do this in release.
#endif
                doGrain(s,
                        x_and_ffts,
                        y,
                        g.second,
                        workData,
                        nForceSteps,
                        yPhase);
                if(g.second == GrainType::FFT) {
                    --nForceSteps;
                }
                ++yPhase;
                s.updatePostGrain(g.second);
                s.grain_counter = 0;
            }
            else {
                s.grain_counter += vectorSz-1;
                break;
            }
        } while (vectorSz);
        
        Assert(nForceSteps==0);
    }
    
    template<template<typename> typename Allocator2, typename WorkData>
    void step(State & s,
              XAndFFTS<T, Allocator2, Tag> const & x_and_ffts,
              Y<T, Tag> & y,
              WorkData * workData) const
    {
        if(unlikely(countPartitions() == 0)) {
            return;
        }
        assert(isValid());
        ++s.grain_counter;
        auto g = nextGrain(s, x_and_ffts.progress);
        assert(g.first >= 0);
        if(unlikely(g.first == 0)) {
            doGrain(s,
                    x_and_ffts,
                    y,
                    g.second,
                    workData,
                    0,
                    0);
            s.updatePostGrain(g.second);
            s.grain_counter = 0;
        }
    }

    void flushToSilence(State & s) const {
        s.flushToSilence();
    }

private:

    std::pair<int, GrainType> nextGrain(State const & s,
                                        unsigned int const x_progress) const {
        
        auto const n_mult_grains_remaining = count_multiplicative_grains - s.grain_number;
        
        if(unlikely(n_mult_grains_remaining < 0)) {
            int distToFFTGrain = partition_size_minus_one - (x_progress & partition_size_minus_one);
            assert(distToFFTGrain >= 0);
            return {distToFFTGrain, GrainType::FFT};
        }
        
        int const dist = getGranularity() - s.grain_counter;
        
        assert(dist >= 0);
        if(unlikely(n_mult_grains_remaining == 0)) {
            return {dist, GrainType::IFFT};
        }
        else {
            Assert(n_mult_grains_remaining > 0);
            return {dist, GrainType::MultiplicationGroup};
        }
    }
    
    template<template<typename> typename Allocator2, typename WorkData>
    void doGrain(State & s,
                 XAndFFTS<T, Allocator2, Tag> const & x_and_ffts,
                 Y<T, Tag> & y,
                 GrainType g,
                 WorkData * workData,
                 int const futureBlock,
                 int const yPhase) const
    {
        if(g == GrainType::FFT) {
            // the fft is implicit (in x_and_ffts)
            //Assert(0 == ((x_and_ffts.progress+1) & partition_size_minus_one)); // make sure 'rythm is good'
        }
        else {
            auto const fft_length = get_fft_length();
            auto const & ffts = find_ffts(x_and_ffts.x_ffts.x_ffts,
                                          fft_length);

            if(g == GrainType::MultiplicationGroup) {
                auto const M = getMultiplicationsGroupMaxSize();
                int offset = M * s.grain_number;
                assert(offset >= 0);
                assert(offset < s.ffts_of_partitionned_h.size());

                int offset_end = std::min(offset + M,
                                          partition_count);
                if(offset == 0) {
                    fft::RealFBins_<Tag, FPT, Allocator>::multiply(s.multiply_add_result.data() /* = */,
                                                                   ffts.get_by_age(futureBlock), /* x */
                                                                   s.ffts_of_partitionned_h[0].data(),
                                                                   1+partition_size_minus_one);
                    offset = 1;
                }
                for(; offset != offset_end; ++offset)
                {
                    fft::RealFBins_<Tag, FPT, Allocator>::multiply_add(s.multiply_add_result.data() /* += */,
                                                                       ffts.get_by_age(futureBlock+offset), /* x  */
                                                                       s.ffts_of_partitionned_h[offset].data(),
                                                                       1+partition_size_minus_one);
                }
            }
            else {
                Assert(g==GrainType::IFFT);

                ffts.fft.inverse(s.multiply_add_result.data(),
                                 workData,
                                 fft_length);

                Assert(s.grain_counter == getGranularity());
                Assert(s.grain_number == count_multiplicative_grains);
                int const blockProgress = (count_multiplicative_grains+1) * getGranularity();
                int const yFutureLocation = y.progress + yPhase + getBlockSize() - blockProgress;
                if constexpr (overlapMode == Overlap::Add) {
                    y.addAssign(yFutureLocation,
                                workData,
                                fft_length);
                }
                else {
                    y.addAssign(yFutureLocation,
                                &workData[fft_length/2],
                                fft_length/2);
                }
            }
        }
    }
    
private:
    int32_t partition_count = 0;
    int32_t count_multiplicative_grains = 0;
    int32_t granularity = 0;
    uint32_t partition_size_minus_one = -1;
    int32_t mult_grp_len = 0; // accessed least frequently
};

}
