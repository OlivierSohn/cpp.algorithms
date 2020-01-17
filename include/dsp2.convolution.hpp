
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
    using Parent::doReset;

    MinSizeRequirement setCoefficients(Algo const & algo, a64::vector<T> coeffs_)
    {
        result.clear();
        result.resize(algo.get_fft_length());
        return doSetCoefficients(algo, std::move(coeffs_));
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
    
    RealSignal result;
};

template <typename Parent>
struct AlgoFFTConvolutionIntermediate : public Parent {
    using State = StateFFTConvolutionIntermediate<typename Parent::State>;
    
    using T = typename Parent::FPT;
    using FPT = T;
    using Tag = typename Parent::FFTTag;
    using Desc = DescFFTConvolutionIntermediate;
    
    using SetupParam = typename Parent::SetupParam;
    
    using RealSignal = typename fft::RealSignal_<Tag, FPT>::type;
    using Algo = typename fft::Algo_<Tag, FPT>;
    using Contexts = fft::Contexts_<Tag, FPT>;
    
    static constexpr auto add_assign = fft::RealSignal_<Tag, FPT>::add_assign;
    
    using Parent::get_fft_length;
    using Parent::compute_convolution;
    using Parent::getBlockSize;
    
    void step(State & s,
              XAndFFTS<T, Tag> const & x_and_ffts,
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
    
    void forceStep(State & s,
                   XAndFFTS<T, Tag> const & x_and_ffts,
                   Y<T, Tag> & y) const
    {
        int const N = getBlockSize();
        
        Assert(y.uProgress+(N-1) < y.y.size());
        add_assign(y.y.begin() + y.uProgress,
                   s.result.begin() + N,
                   N);
        
        auto const fft_length = get_fft_length();

        auto const & ffts = x_and_ffts.find_ffts(fft_length);
        
        auto const & frequencies = compute_convolution(s, ffts);
        
        ffts.fft.inverse(frequencies, s.result, fft_length);
        
        add_assign(y.y.begin() + y.uProgress,
                   s.result.begin(),
                   N);
    }
    
    void flushToSilence(State & s) const {
        s.flushToSilence();
    }
};

template<typename Parent>
struct corresponding_legacy_dsp<AlgoFFTConvolutionIntermediate<Parent>> {
    using type = FFTConvolutionIntermediate<corresponding_legacy_dsp_t<Parent>>;
};


/*
 */


template <typename T, typename Tag>
struct AlgoFFTConvolutionCRTP;

template <typename T, typename Tag>
struct StateFFTConvolutionCRTP {
    using Algo = AlgoFFTConvolutionCRTP<T, Tag>;

    using FPT = T;
    using FFTTag = Tag;
    
    using RealSignal = typename fft::RealSignal_<Tag, FPT>::type;
    static constexpr auto makeRealSignal = fft::RealSignal_<Tag, FPT>::make;
    
    static constexpr auto scale = fft::RealFBins_<Tag, FPT>::scale;
    
    using CplxFreqs = typename fft::RealFBins_<Tag, FPT>::type;
        
    MinSizeRequirement doSetCoefficients(Algo const & algo, a64::vector<T> coeffs_)
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
            fft.forward(coeffs.begin(), fft_of_h, fft_length);
            
            auto factor = scaleFactor<FFTAlgo>(static_cast<FPT>(algo.get_fft_length()));
            scale(fft_of_h, factor);
            
            Assert(fft_length>1);
        }
        return {
            0, // x block size
            static_cast<int>(fft_length/2), // y block size
            0, // no y anticipated writes
            {
                {fft_length, 1}
            }
        };
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

template <typename T, typename Tag>
struct AlgoFFTConvolutionCRTP {

    using State = StateFFTConvolutionCRTP<T, Tag>;
    
    using FPT = T;
    using FFTTag = Tag;
    
    using CplxFreqs = typename fft::RealFBins_<Tag, FPT>::type;
    static constexpr auto multiply = fft::RealFBins_<Tag, FPT>::multiply;
    
    using FFTAlgo = typename fft::Algo_<Tag, FPT>;
    
    using SetupParam = FFTConvolutionCRTPSetupParam;
    
    void setup(SetupParam const & p) {
        N = p.blockSize;
        if(N>0 && !is_power_of_two(N)) {
            throw std::runtime_error("non power of 2");
        }
        auto fft_length = 2*N;
        work.resize(fft_length);
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
    
    auto const & compute_convolution(State & s,
                                     typename XAndFFTS<T, Tag>::FFTs const & ffts) const
    {
        auto const & fft_of_x = ffts.get_backward(0);
        
        multiply(work, fft_of_x, s.fft_of_h);
        
        return work;
    }
    
private:
    int N;
    mutable CplxFreqs work;
};

template <typename T, typename Tag>
struct corresponding_legacy_dsp<AlgoFFTConvolutionCRTP<T, Tag>> {
    using type = FFTConvolutionCRTP<T, Tag>;
};


/*
 */


template <typename T, typename Tag>
struct AlgoPartitionnedFFTConvolutionCRTP;

template <typename T, typename Tag>
struct StatePartitionnedFFTConvolutionCRTP {
    using FPT = T;
    using FFTTag = Tag;
    using Algo = AlgoPartitionnedFFTConvolutionCRTP<T, Tag>;
    
    using FFTAlgo = typename fft::Algo_<Tag, FPT>;
    
    using RealSignal = typename fft::RealSignal_<Tag, FPT>::type;
    using Signal_value_type = typename fft::RealSignal_<Tag, FPT>::value_type;
    
    using CplxFreqs = typename fft::RealFBins_<Tag, FPT>::type;
    
    static constexpr auto scale = fft::RealFBins_<Tag, FPT>::scale;
    
    MinSizeRequirement doSetCoefficients(Algo const & algo,
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
                using FFTAlgo = typename fft::Algo_<Tag, FPT>;
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
                    fft.forward(coeffs_slice.begin(), fft_of_partitionned_h, fft_length);
                    scale(fft_of_partitionned_h, factor);
                }
            }
            Assert(it_coeffs == coeffs_.end());
            
            Assert(fft_length > 1);
        }
        return {
            0, // x block size
            static_cast<int>(fft_length/2), // y block size
            0, // no y anticipated writes
            {
                {fft_length, n_partitions}
            }
        };
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

template <typename T, typename Tag>
struct AlgoPartitionnedFFTConvolutionCRTP {

    using SetupParam = PartitionnedFFTConvolutionCRTPSetupParam;
    using FPT = T;
    using FFTTag = Tag;
    using State = StatePartitionnedFFTConvolutionCRTP<T, Tag>;
    
    using RealSignal = typename fft::RealSignal_<Tag, FPT>::type;
    using Signal_value_type = typename fft::RealSignal_<Tag, FPT>::value_type;
    
    using CplxFreqs = typename fft::RealFBins_<Tag, FPT>::type;
    static constexpr auto multiply = fft::RealFBins_<Tag, FPT>::multiply;
    static constexpr auto multiply_add = fft::RealFBins_<Tag, FPT>::multiply_add;
    
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

        work.resize(get_fft_length());
    }
    
    auto countPartitions() const { return partition_count; }
    
protected:
    auto const & compute_convolution(State & s,
                                     typename XAndFFTS<T, Tag>::FFTs const & ffts) const
    {
        int index = 0;
        
        struct F {
            void operator()(CplxFreqs const & fft_of_delayed_x) {
                auto const & fft_of_partitionned_h = ffts_of_partitionned_h[index];
                if(0 == index) {
                    multiply(work /* = */, fft_of_delayed_x, /* * */ fft_of_partitionned_h);
                }
                else {
                    multiply_add(work /* += */, fft_of_delayed_x, /* * */ fft_of_partitionned_h);
                }
                ++index;
            }
            int & index;
            std::vector<CplxFreqs> & ffts_of_partitionned_h;
            CplxFreqs & work;
        } f{index, s.ffts_of_partitionned_h, work};
        
        ffts.for_some_bkwd(s.ffts_of_partitionned_h.size(), f);
        
        return work;
    }
    
private:
    int partition_size = -1;
    int partition_count = 0;
    mutable CplxFreqs work;
};

template <typename T, typename Tag>
struct corresponding_legacy_dsp<AlgoPartitionnedFFTConvolutionCRTP<T, Tag>> {
    using type = PartitionnedFFTConvolutionCRTP<T, Tag>;
};

}
