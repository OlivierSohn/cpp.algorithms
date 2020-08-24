namespace imajuscule::audio {

/*
 Notations for complexity:
 
 H : length of impulse response
 
 A : number of frames computed during one audio callback
 
 PART_N : Number of partitions
 PART_L : Length of one partition
 ( H == PART_N * PART_L )
 
 The computations are based on the fact that an fft of an S-long signal costs S*lg(S)
 */

template<typename Algo>
auto scaleFactor(typename Algo::FPT fft_length) -> typename Algo::FPT {
    return 1. / (Algo::scale * Algo::scale * fft_length);
}

struct DescFFTConvolutionIntermediate {
    static constexpr int nCoefficientsFadeIn = 0;
    static constexpr bool has_subsampling = false;
    static constexpr bool step_can_error = false;
};

struct StepsData {
    StepsData(int nBlockSteps,
              int yPhase)
    : nBlockSteps(nBlockSteps)
    , yPhase(yPhase)
    {}
    
    int nBlockSteps;
    int yPhase;
};

template <typename Parent>
struct AlgoFFTConvolutionIntermediate : public Parent {
    using State = typename Parent::State;
    
    using T = typename Parent::FPT;
    using FPT = T;
    using Tag = typename Parent::Tag;
    
    using FFTAlgo = fft::Algo_<Tag, T>;
    
    template<typename TT>
    using Allocator = typename Parent::template Allocator<TT>;
    
    using Desc = DescFFTConvolutionIntermediate;
    
    using SetupParam = typename Parent::SetupParam;
    
    using RealSignal = typename fft::RealSignal_<Tag, FPT>::type;
    using Contexts = fft::Contexts_<Tag, FPT>;
    
    using Parent::compute_convolution;
    using Parent::doDephaseStep;
    using Parent::get_fft_length;
    using Parent::getBlockSize;

    void dephaseStep(State & s, int x_progress) const {
        doDephaseStep(s, x_progress);
    }
    
    StepsData getStepsData(int const x_progress,
                           int const vectorSz) const
    {
        auto const N = getBlockSize();
        int const h = static_cast<unsigned int>(x_progress + 1) & ~(N-1);
        int const l = static_cast<unsigned int>(x_progress + 1 - vectorSz) & ~(N-1);
        int const log2Blocksize = count_trailing_zeroes(N);
        
        return {
            (h-l) >> log2Blocksize,
            static_cast<int>((N - ((static_cast<unsigned int>(x_progress + 2 - vectorSz)) & (N-1))) & (N-1))
        };
        
        /*
        int nForceSteps2 = x_and_ffts.countFftsSinceForhalfSize(vectorSz,
                                                               N);
        Assert(nForceSteps==nForceSteps2);
        */
    }
    
    template<template<typename> typename Allocator2, typename WorkData>
    void stepVectorized(State & s,
                        XAndFFTS<T, Allocator2, Tag> const & x_and_ffts,
                        Y<T, Tag> & y,
                        WorkData * workData,
                        int const vectorSz) const
    {
        if(unlikely(s.isZero())) {
            return;
        }
        auto stepsData = getStepsData(x_and_ffts.progress,
                                      vectorSz);
        if(!stepsData.nBlockSteps) {
            return;
        }
        forceSteps(s,
                   x_and_ffts,
                   y,
                   workData,
                   stepsData);
    }
    
    template<template<typename> typename Allocator2, typename WorkData>
    void step(State & s,
              XAndFFTS<T, Allocator2, Tag> const & x_and_ffts,
              Y<T, Tag> & y,
              WorkData * workData) const
    {
        if(s.isZero()) {
            return;
        }

        unsigned int const N = getBlockSize();
        Assert(!N || is_power_of_two(N));
        
        if((static_cast<unsigned int>(x_and_ffts.progress + 1) & (N-1)) == 0) {
            forceSteps(s,
                       x_and_ffts,
                       y,
                       workData,
                       {1, 0});
        }
    }
    
    template<template<typename> typename Allocator2, typename WorkData>
    void forceSteps(State const & s,
                    XAndFFTS<T, Allocator2, Tag> const & x_and_ffts,
                    Y<T, Tag> & y,
                    WorkData * workData,
                    StepsData const & d) const
    {
        Assert(d.nBlockSteps);
        
        auto const fft_length = get_fft_length();

        auto const & ffts = find_ffts(x_and_ffts.x_ffts.x_ffts,
                                      fft_length);
        
        for(int i=0; i<d.nBlockSteps; ++i) {
            compute_convolution(s, ffts, workData, d.nBlockSteps-i-1);
            
            if constexpr (FFTAlgo::inplace_dft)
            {
                ffts.fft.inverse(&workData[0],
                                 fft_length);

                if(i==0 && d.yPhase == 0) {
                    if constexpr (overlapMode == Overlap::Add) {
                        y.inplace_addAssignPresent(&workData[0],
                                                   0,
                                                   fft_length);
                    }
                    else {
                        y.inplace_addAssignPresent(&workData[0],
                                                   fft_length/2,
                                                   fft_length/2);
                    }
                }
                else {
                    if constexpr (overlapMode == Overlap::Add) {
                        y.inplace_addAssign(y.progress + d.yPhase + i*getBlockSize(),
                                            &workData[0],
                                            0,
                                            fft_length);
                    }
                    else {
                        y.inplace_addAssign(y.progress + d.yPhase + i*getBlockSize(),
                                            &workData[0],
                                            fft_length/2,
                                            fft_length/2);
                    }
                }
            }
            else {
                ffts.fft.inverse(&workData[0],
                                 &workData[fft_length],
                                 fft_length);
                
                if(i==0 && d.yPhase == 0) {
                    if constexpr (overlapMode == Overlap::Add) {
                        y.addAssignPresent(&workData[fft_length],
                                           fft_length);
                    }
                    else {
                        y.addAssignPresent(&workData[fft_length+fft_length/2],
                                           fft_length/2);
                    }
                }
                else {
                    if constexpr (overlapMode == Overlap::Add) {
                        y.addAssign(y.progress + d.yPhase + i*getBlockSize(),
                                    &workData[fft_length],
                                    fft_length);
                    }
                    else {
                        y.addAssign(y.progress + d.yPhase + i*getBlockSize(),
                                    &workData[fft_length+fft_length/2],
                                    fft_length/2);
                    }
                }
            }
        }
    }
    
    void flushToSilence(State & s) const {
        s.flushToSilence();
    }
};


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
template <typename T, template<typename> typename Allocator, typename Tag>
struct AlgoFFTConvolutionCRTP;

template <typename T, template<typename> typename Allocator, typename FFTTag>
struct StateFFTConvolutionCRTP {
    using Algo = AlgoFFTConvolutionIntermediate<AlgoFFTConvolutionCRTP<T, Allocator, FFTTag>>;

    using Desc = DescFFTConvolutionIntermediate;

    using FPT = T;
    using Tag = FFTTag;
    
    using RealSignal = typename fft::RealSignal_<Tag, FPT>::type;
    static constexpr auto makeRealSignal = fft::RealSignal_<Tag, FPT>::make;
    
    static constexpr auto scale = fft::RealFBins_<Tag, FPT, Allocator>::scale;
    
    using CplxFreqs = typename fft::RealFBins_<Tag, FPT, Allocator>::type;

    static int getAllocationSz_SetCoefficients(typename Algo::SetupParam const & p) {
        auto fft_length = 2 * p.blockSize;
        return fft_length;
    }
    void setCoefficients(Algo const & algo, a64::vector<T> coeffs_)
    {
        auto fft_length = algo.get_fft_length();
        if(fft_length) {
            fft_of_h.resize(fft_length);
            
            // pad impulse response with 0
            
            coeffs_.resize(fft_length, {});
            
            // compute fft of padded impulse response
            auto coeffs = makeRealSignal(std::move(coeffs_));
            
            using FFTAlgo = typename fft::Algo_<Tag, FPT>;
            using Contexts = fft::Contexts_<Tag, FPT>;
            FFTAlgo fft(Contexts::getInstance().getBySize(fft_length));
            fft.forward(coeffs.begin(), fft_of_h.data(), fft_length);
            
            auto factor = scaleFactor<FFTAlgo>(static_cast<FPT>(algo.get_fft_length()));
            scale(fft_of_h, factor);
            
            Assert(fft_length>1);
        }
    }

    void flushToSilence() {
    }

    template<typename F>
    void onContextFronteer(F f) {
    }
    
    bool isZero() const {
        return fft_of_h.empty();
    }
    
    void reset() {
        fft_of_h.clear();
    }

    double getEpsilon(Algo const & algo) const {
        return fft::getFFTEpsilon<FPT>(algo.get_fft_length()) + 2 * std::numeric_limits<FPT>::epsilon();
    }
    void logComputeState(Algo const & algo, std::ostream & os) const {
        os << "Ffts, 1 block of size " << algo.getBlockSize() << std::endl;
    }
    
public:
    CplxFreqs fft_of_h;
};

template <typename T, template<typename> typename Alloc, typename FFTTag>
struct AlgoFFTConvolutionCRTP {
    template<typename TT>
    using Allocator = Alloc<TT>;
    
    using State = StateFFTConvolutionCRTP<T, Alloc, FFTTag>;
    
    using FPT = T;
    using Tag = FFTTag;
    
    using CplxFreqsValueType = typename fft::RealFBins_<Tag, FPT, Allocator>::type::value_type;
    
    using SetupParam = FFTConvolutionCRTPSetupParam;
    
    void setup(SetupParam const & p) {
        N = p.blockSize;
        if(N>0 && !is_power_of_two(N)) {
            throw std::runtime_error("non power of 2");
        }
    }
    
    void doDephaseStep(State & s,
                       int x_progress) const {
        // dephasing occurs right after initialization when all buffers are 0s.
        // since there is no other state that work buffers (that are 0s) there is no need to do anything here.
    }

    bool isValid() const { return true; }
    bool handlesCoefficients() const { return N > 0; }

    auto get_fft_length() const { return 2 * N; }
    auto getBlockSize() const { return N; }
    Latency getLatency() const {
        Assert(handlesCoefficients());
        return Latency(N-1);
    }
    auto countPartitions() const { return 1; }
    
    static auto get_fft_length(int n) {
        auto N_nonzero_y = 2 * n;
        return ceil_power_of_two(N_nonzero_y);
    }
    
protected:
    template<template<typename> typename Allocator2>
    void compute_convolution(State const & s,
                             FFTs<T, Allocator2, Tag> const & ffts,
                             CplxFreqsValueType * workData,
                             int const futureBlock) const
    {
        fft::RealFBins_<Tag, FPT, Allocator2>::multiply(workData,
                                                        ffts.get_by_age(futureBlock),
                                                        s.fft_of_h.data(),
                                                        N);
    }
    
private:
    int N;
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
 *
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
template <typename T, template<typename> typename Allocator, typename Tag>
struct AlgoPartitionnedFFTConvolutionCRTP;

template <typename T, template<typename> typename Alloc, typename FFTTag>
struct StatePartitionnedFFTConvolutionCRTP {

    template<typename TT>
    using Allocator = Alloc<TT>;

    using Desc = DescFFTConvolutionIntermediate;

    using FPT = T;
    using Tag = FFTTag;
    using Algo = AlgoFFTConvolutionIntermediate<AlgoPartitionnedFFTConvolutionCRTP<T, Alloc, Tag>>;
    
    using FFTAlgo = typename fft::Algo_<Tag, FPT>;
    
    using RealSignal = typename fft::RealSignal_<Tag, FPT>::type;
    using Signal_value_type = typename fft::RealSignal_<Tag, FPT>::value_type;
    
    using CplxFreqs = typename fft::RealFBins_<Tag, FPT, Allocator>::type;
    
    static constexpr auto scale = fft::RealFBins_<Tag, FPT, Allocator>::scale;
    
    static int getAllocationSz_SetCoefficients(typename Algo::SetupParam const & p) {
        auto const fft_length = 2 * p.partition_size;
        return fft_length * p.partition_count;
    }
    void setCoefficients(Algo const & algo,
                           a64::vector<T> coeffs_)
    {
        int const partition_size = algo.getBlockSize();
        auto const n_partitions = countPartitions(coeffs_.size(),
                                                  partition_size);
        if(n_partitions != algo.countPartitions()) {
            throw std::logic_error("inconsistent count of partitions");
        }
        // if one partition is partial, pad it with zeros
        coeffs_.resize(n_partitions * partition_size,
                       {});

        ffts_of_partitionned_h.resize(n_partitions);

        auto const fft_length = algo.get_fft_length();
        if(fft_length) {
            for(auto & fft_of_partitionned_h : ffts_of_partitionned_h) {
                fft_of_partitionned_h.resize(fft_length);
            }
            
            // compute fft of padded impulse response
            
            auto it_coeffs = coeffs_.begin();
            {
                using Contexts = fft::Contexts_<Tag, FPT>;
                FFTAlgo fft(Contexts::getInstance().getBySize(fft_length));
                
                auto const factor = scaleFactor<FFTAlgo>(static_cast<FPT>(fft_length));
                RealSignal coeffs_slice(fft_length, Signal_value_type(0)); // initialize with zeros (second half is padding)
                for(auto & fft_of_partitionned_h : ffts_of_partitionned_h) {
                    auto end_coeffs = it_coeffs + algo.getBlockSize();
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
    }

    void flushToSilence() {
    }

    template<typename F>
    void onContextFronteer(F f) {
    }

    void reset() {
        ffts_of_partitionned_h.clear();
    }

    bool isZero() const {
        return ffts_of_partitionned_h.empty();
    }

    double getEpsilon(Algo const & algo) const {
        return algo.countPartitions() * (fft::getFFTEpsilon<FPT>(algo.get_fft_length()) + 2 * std::numeric_limits<FPT>::epsilon());
    }
    void logComputeState(Algo const & algo, std::ostream & os) const {
        os << "FFTs, " << algo.countPartitions() << " partitions of sizes " << algo.getBlockSize() << " each" << std::endl;
    }
    
    std::vector<CplxFreqs> ffts_of_partitionned_h;
};

template <typename T, template<typename> typename Alloc, typename FFTTag>
struct AlgoPartitionnedFFTConvolutionCRTP {
    template<typename TT>
    using Allocator = Alloc<TT>;
    
    using Desc = DescFFTConvolutionIntermediate;
    using SetupParam = PartitionnedFFTConvolutionCRTPSetupParam;
    using FPT = T;
    using Tag = FFTTag;
    using State = StatePartitionnedFFTConvolutionCRTP<T, Alloc, Tag>;
    
    using RealSignal = typename fft::RealSignal_<Tag, FPT>::type;
    using Signal_value_type = typename fft::RealSignal_<Tag, FPT>::value_type;
    
    using CplxFreqsValueType = typename fft::RealFBins_<Tag, FPT, Allocator>::type::value_type;
    using CplxFreqs = typename fft::RealFBins_<Tag, FPT, Allocator>::type;
    
    auto get_fft_length() const { assert(partition_size > 0); return 2 * partition_size; }
    auto get_fft_length(int) const { return get_fft_length(); }
    auto getBlockSize() const { return partition_size; }
    Latency getLatency() const {
        Assert(handlesCoefficients());
        return Latency(partition_size-1);
    }
    
    bool isValid() const { return true; }
    bool handlesCoefficients() const { return partition_size > 0; }

    void setup(SetupParam const & p) {
        partition_size = p.partition_size;
        partition_count = p.partition_count;
        
        assert(partition_size > 0);
        assert(is_power_of_two(partition_size));
    }

    void doDephaseStep(State & s,
                       int x_progress) const {
        // dephasing occurs right after initialization when all buffers are 0s.
        // since there is no other state that work buffers (that are 0s) there is no need to do anything here.
    }

    auto countPartitions() const { return partition_count; }
    
protected:
    
    template<template<typename> typename Allocator2>
    void compute_convolution(State const & s,
                             FFTs<T, Allocator2, Tag> const & ffts,
                             CplxFreqsValueType * workData,
                             int futureBlock) const
    {
        int index = 0;
        using CplxFreqs2 = typename fft::RealFBins_<Tag, FPT, Allocator2>::type;
        struct F {
            void operator()(CplxFreqs2 const & fft_of_delayed_x) {
                auto const & fft_of_partitionned_h = ffts_of_partitionned_h[index];
                if(0 == index) {
                    fft::RealFBins_<Tag, FPT, Allocator>::multiply(workData /* = */,
                                                                   fft_of_delayed_x.data(), /* * */ fft_of_partitionned_h.data(),
                                                                   partition_size);
                }
                else {
                    fft::RealFBins_<Tag, FPT, Allocator>::multiply_add(workData /* += */,
                                                                       fft_of_delayed_x.data(), /* * */ fft_of_partitionned_h.data(),
                                                                       partition_size);
                }
                ++index;
            }
            int & index;
            std::vector<typename State::CplxFreqs> const & ffts_of_partitionned_h;
            CplxFreqsValueType * workData;
            int partition_size;
        } f{index, s.ffts_of_partitionned_h, workData, partition_size};
        
        ffts.for_some_recent(futureBlock, s.ffts_of_partitionned_h.size(), f);
    }
    
private:
    int partition_size = -1;
    int partition_count = 0;
};

}
