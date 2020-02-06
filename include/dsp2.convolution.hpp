
namespace imajuscule {

template<typename Algo>
auto scaleFactor(typename Algo::FPT fft_length) -> typename Algo::FPT {
    return 1. / (Algo::scale * Algo::scale * fft_length);
}

struct DescFFTConvolutionIntermediate {
    static constexpr int nCoefficientsFadeIn = 0;
    static constexpr bool has_subsampling = false;
    static constexpr bool step_can_error = false;
};

template <typename Parent>
struct AlgoFFTConvolutionIntermediate;

template <typename Parent>
struct StateFFTConvolutionIntermediate : public Parent {
    
    using FPT = typename Parent::FPT;
    using T = FPT;
    using Tag = typename Parent::FFTTag;
    using Desc = DescFFTConvolutionIntermediate;
    using Algo = AlgoFFTConvolutionIntermediate<typename Parent::Algo>;
    
    using RealSignal = typename fft::RealSignal_<Tag, FPT>::type;

    using Parent::doSetCoefficients;
    using Parent::doFlushToSilence;
    using Parent::doOnContextFronteer;
    using Parent::doReset;
    using Parent::doGetAllocationSz_SetCoefficients;

    static int getAllocationSz_SetCoefficients(typename Algo::SetupParam const & p) {
        return doGetAllocationSz_SetCoefficients(p);
    }
    void setCoefficients(Algo const & algo, a64::vector<T> coeffs_)
    {
        result.clear();
        result.resize(algo.get_fft_length());
        doSetCoefficients(algo, std::move(coeffs_));
    }
    
    template<typename F>
    void onContextFronteer(F f) {
        doOnContextFronteer(f);
    }
    
    void reset() {
        result.clear();
        doReset();
    }
    
    void flushToSilence() {
        Parent::doFlushToSilence();

        fft::RealSignal_<Tag, FPT>::zero(result);
    }

    bool isZero() const {
        return result.empty();
    }
    
    // mutable because used to store results of previous step, to apply them in the next step
    mutable RealSignal result;
};

template <typename Parent>
struct AlgoFFTConvolutionIntermediate : public Parent {
    using State = StateFFTConvolutionIntermediate<typename Parent::State>;
    
    using T = typename Parent::FPT;
    using FPT = T;
    using Tag = typename Parent::FFTTag;
    
    template<typename TT>
    using Allocator = typename Parent::template Allocator<TT>;
    
    using Desc = DescFFTConvolutionIntermediate;
    
    using SetupParam = typename Parent::SetupParam;
    
    using RealSignal = typename fft::RealSignal_<Tag, FPT>::type;
    using Contexts = fft::Contexts_<Tag, FPT>;
    
    static constexpr auto add_assign = fft::RealSignal_<Tag, FPT>::add_assign;
    
    using Parent::compute_convolution;
    using Parent::doDephaseSteps;
    using Parent::get_fft_length;
    using Parent::getBlockSize;

    void dephaseSteps(State & s, int n_steps) const {
        doDephaseSteps(s, n_steps);
    }

    template<template<typename> typename Allocator2, typename WorkCplxFreqs>
    void step(State & s,
              XAndFFTS<T, Allocator2, Tag, WorkCplxFreqs> const & x_and_ffts,
              Y<T, Tag> & y) const
    {
        if(s.isZero()) {
            return;
        }

        auto const N = getBlockSize();
        if(unlikely(!N)) {
            return;
        }
        auto fft_length = get_fft_length();
        Assert(fft_length == 2*N);
        Assert(is_power_of_two(N));
        
        if((x_and_ffts.fftsHalfSizesBitsToCompute & static_cast<uint32_t>(N)) != N) {
            return;
        }
        
        Assert((static_cast<uint32_t>(x_and_ffts.progress) & static_cast<uint32_t>(getBlockSize()-1)) == 0);
        Assert(s.result.size() == 2*N);
        
        forceStep(s,
                  x_and_ffts,
                  y);
    }
    
    template<template<typename> typename Allocator2, typename WorkCplxFreqs>
    void forceStep(State const & s,
              XAndFFTS<T, Allocator2, Tag, WorkCplxFreqs> const & x_and_ffts,
              Y<T, Tag> & y) const
    {
        int const N = getBlockSize();
        
        Assert(y.uProgress+(N-1) < y.y.size());
        add_assign(y.y.begin() + y.uProgress,
                   s.result.begin() + N,
                   N);
        
        auto const fft_length = get_fft_length();

        auto const & ffts = x_and_ffts.find_ffts(fft_length);
        
        auto * workData = x_and_ffts.work.data();
        compute_convolution(s, ffts, workData);
        
        ffts.fft.inverse(workData, s.result, fft_length);
        
        add_assign(y.y.begin() + y.uProgress,
                   s.result.begin(),
                   N);
    }
    
    void flushToSilence(State & s) const {
        s.flushToSilence();
    }
};


/*
 */


template <typename T, template<typename> typename Allocator, typename Tag>
struct AlgoFFTConvolutionCRTP;

template <typename T, template<typename> typename Allocator, typename Tag>
struct StateFFTConvolutionCRTP {
    using Algo = AlgoFFTConvolutionCRTP<T, Allocator, Tag>;

    using FPT = T;
    using FFTTag = Tag;
    
    using RealSignal = typename fft::RealSignal_<Tag, FPT>::type;
    static constexpr auto makeRealSignal = fft::RealSignal_<Tag, FPT>::make;
    
    static constexpr auto scale = fft::RealFBins_<Tag, FPT, Allocator>::scale;
    
    using CplxFreqs = typename fft::RealFBins_<Tag, FPT, Allocator>::type;

    static int doGetAllocationSz_SetCoefficients(typename Algo::SetupParam const & p) {
        auto fft_length = 2 * p.blockSize;
        return fft_length;
    }
    void doSetCoefficients(Algo const & algo, a64::vector<T> coeffs_)
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
    
    template<typename F>
    void doOnContextFronteer(F f) {
    }
    
    void doReset() {
        fft_of_h.clear();
    }

    double getEpsilon(Algo const & algo) const {
        return fft::getFFTEpsilon<FPT>(algo.get_fft_length()) + 2 * std::numeric_limits<FPT>::epsilon();
    }
    void logComputeState(Algo const & algo, std::ostream & os) const {
        os << "Ffts, 1 block of size " << algo.getBlockSize() << std::endl;
    }
protected:
    void doFlushToSilence() {
        // nothing to do : no member contains state related to a past signal, except
        // for fft_of_x but it will be overwritten next time we use it in compute_convolution
    }
    
public:
    CplxFreqs fft_of_h;
};

template <typename T, template<typename> typename Alloc, typename Tag>
struct AlgoFFTConvolutionCRTP {
    template<typename TT>
    using Allocator = Alloc<TT>;
    
    using State = StateFFTConvolutionCRTP<T, Alloc, Tag>;
    
    using FPT = T;
    using FFTTag = Tag;
    
    using CplxFreqsValueType = typename fft::RealFBins_<Tag, FPT, Allocator>::type::value_type;
    
    using FFTAlgo = typename fft::Algo_<Tag, FPT>;
    
    using SetupParam = FFTConvolutionCRTPSetupParam;
    
    void setup(SetupParam const & p) {
        N = p.blockSize;
        if(N>0 && !is_power_of_two(N)) {
            throw std::runtime_error("non power of 2");
        }
    }
    

    void doDephaseSteps(State & s,
                        int n_steps) const {
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
                             CplxFreqsValueType * workData) const
    {
        fft::RealFBins_<Tag, FPT, Allocator2>::multiply(workData,
                                                        ffts.get_by_age(0),
                                                        s.fft_of_h.data(),
                                                        N);
    }
    
private:
    int N;
};

/*
 */


template <typename T, template<typename> typename Allocator, typename Tag>
struct AlgoPartitionnedFFTConvolutionCRTP;

template <typename T, template<typename> typename Alloc, typename Tag>
struct StatePartitionnedFFTConvolutionCRTP {

    template<typename TT>
    using Allocator = Alloc<TT>;

    using FPT = T;
    using FFTTag = Tag;
    using Algo = AlgoPartitionnedFFTConvolutionCRTP<T, Alloc, Tag>;
    
    using FFTAlgo = typename fft::Algo_<Tag, FPT>;
    
    using RealSignal = typename fft::RealSignal_<Tag, FPT>::type;
    using Signal_value_type = typename fft::RealSignal_<Tag, FPT>::value_type;
    
    using CplxFreqs = typename fft::RealFBins_<Tag, FPT, Allocator>::type;
    
    static constexpr auto scale = fft::RealFBins_<Tag, FPT, Allocator>::scale;
    
    static int doGetAllocationSz_SetCoefficients(typename Algo::SetupParam const & p) {
        auto const fft_length = 2 * p.partition_size;
        return fft_length * p.partition_count;
    }
    void doSetCoefficients(Algo const & algo,
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

    template<typename F>
    void doOnContextFronteer(F f) {
    }

    void doReset() {
        ffts_of_partitionned_h.clear();
    }
    
    double getEpsilon(Algo const & algo) const {
        return algo.countPartitions() * (fft::getFFTEpsilon<FPT>(algo.get_fft_length()) + 2 * std::numeric_limits<FPT>::epsilon());
    }
    void logComputeState(Algo const & algo, std::ostream & os) const {
        os << "FFTs, " << algo.countPartitions() << " partitions of sizes " << algo.getBlockSize() << " each" << std::endl;
    }

protected:
    void doFlushToSilence() {
    }
    
public:
    std::vector<CplxFreqs> ffts_of_partitionned_h;
};

template <typename T, template<typename> typename Alloc, typename Tag>
struct AlgoPartitionnedFFTConvolutionCRTP {
    template<typename TT>
    using Allocator = Alloc<TT>;
    
    using SetupParam = PartitionnedFFTConvolutionCRTPSetupParam;
    using FPT = T;
    using FFTTag = Tag;
    using State = StatePartitionnedFFTConvolutionCRTP<T, Alloc, Tag>;
    
    using RealSignal = typename fft::RealSignal_<Tag, FPT>::type;
    using Signal_value_type = typename fft::RealSignal_<Tag, FPT>::value_type;
    
    using CplxFreqsValueType = typename fft::RealFBins_<Tag, FPT, Allocator>::type::value_type;
    using CplxFreqs = typename fft::RealFBins_<Tag, FPT, Allocator>::type;
    
    using FFTAlgo = typename fft::Algo_<Tag, FPT>;
    
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

    void doDephaseSteps(State & s,
                        int n_steps) const {
        // dephasing occurs right after initialization when all buffers are 0s.
        // since there is no other state that work buffers (that are 0s) there is no need to do anything here.
    }

    auto countPartitions() const { return partition_count; }
    
protected:
    
    template<template<typename> typename Allocator2>
    void compute_convolution(State const & s,
                             FFTs<T, Allocator2, Tag> const & ffts,
                             CplxFreqsValueType * workData) const
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
        
        ffts.for_some_recent(s.ffts_of_partitionned_h.size(), f);
    }
    
private:
    int partition_size = -1;
    int partition_count = 0;
};

}
