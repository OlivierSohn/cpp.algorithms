
namespace imajuscule {

struct ConvReverbOptimizationReportSuccess {
    std::string optimizationReport;
};

using ConvReverbOptimizationReport = Either<std::string, ConvReverbOptimizationReportSuccess>;

enum class ReverbType {
    Offline,
    Realtime_Synchronous,
    //Realtime_Synchronous_Subsampled,
    Realtime_Asynchronous_Legacy,
    Realtime_Asynchronous
};

bool constexpr isAsync(ReverbType r) {
    return
    (r == ReverbType::Realtime_Asynchronous_Legacy) ||
    (r == ReverbType::Realtime_Asynchronous);
}

static inline std::string toJustifiedString(ReverbType t) {
    switch(t) {
        case ReverbType::Offline :
            return "Offline                        ";
        case ReverbType::Realtime_Synchronous :
            return "Realtime_Synchronous           ";
            /*case ReverbType::Realtime_Synchronous_Subsampled :
             return "Realtime_Synchronous_Subsampled";
             */
        case ReverbType::Realtime_Asynchronous_Legacy :
            return "Realtime_Asynchronous_Legacy   ";
        case ReverbType::Realtime_Asynchronous :
            return "Realtime_Asynchronous          ";
    }
    return "?";
}

template<typename C>
int countScales(C & rev) {
    auto & lateHandler = rev.getB();
    {
        auto & inner = lateHandler.getB().getInner().getInner();
        if(inner.isZero()) {
            return 1;
        }
    }
    {
        auto & inner = lateHandler.getB().getInner().getInner().getB().getInner().getInner();
        if(inner.isZero()) {
            return 2;
        }
    }
    {
        auto & inner = lateHandler.getB().getInner().getInner().getB().getInner().getInner().getB().getInner().getInner();
        if(inner.isZero()) {
            return 3;
        }
    }
    static_assert(4==nMaxScales);
    return 4;
}

template<int nAudioOut, ReverbType reverbType, PolicyOnWorkerTooSlow OnWorkerTooSlow>
struct Reverbs {
    static constexpr auto nOut = nAudioOut;
    static_assert(nAudioOut > 0);
    static constexpr auto nEars = nAudioOut;
    
    using Tag = fft::Fastest;
    
    template<typename TT>
    using Allocator = monotonic::aP::Alloc<TT>;
    
    using FPT = double;
    using T = FPT;

    using FBinsAllocator = typename fft::RealFBins_<Tag, FPT, Allocator>::type::allocator_type;
    using MemResource = MemResource<FBinsAllocator>;
    
    using State =
    std::conditional_t< reverbType==ReverbType::Offline,
    StateOptimizedFIRFilter<FPT, Allocator, Tag>,
    
    std::conditional_t< reverbType==ReverbType::Realtime_Synchronous,
    StateZeroLatencyScaledFineGrainedPartitionnedConvolution<FPT, Allocator, Tag>,
    
    //std::conditional_t< reverbType==ReverbType::Realtime_Synchronous_Subsampled,
    //  ZeroLatencyScaledFineGrainedPartitionnedSubsampledConvolution<double, Allocator, Tag>,
    
    std::conditional_t< reverbType==ReverbType::Realtime_Asynchronous_Legacy,
    StateZeroLatencyScaledAsyncConvolution<FPT, Allocator, Tag, OnWorkerTooSlow>,
    
    std::conditional_t< reverbType==ReverbType::Realtime_Asynchronous,
    StateZeroLatencyScaledAsyncConvolutionOptimized<FPT, Allocator, Tag, OnWorkerTooSlow>,
    
    void
    
    >>>>/*>*/;
    
    using SetupParam = typename State::SetupParam;
    using PartitionAlgo = PartitionAlgo<SetupParam, FPT, Tag>;
    using WorkCplxFreqs = typename State::WorkCplxFreqs;
    using Algo = typename State::Algo;

    static int getAllocationSz(SetupParam const & p,
                               int const n_channels) {
        return n_channels * (State::getAllocationSz_Setup(p) + State::getAllocationSz_SetCoefficients(p));
    }

    // for wir files of wave, it seems the order is by "ears" then by "source":
    // ear Left source 1
    // ...
    // ear Left source N
    // ear Right source 1
    // ...
    // ear Right source N
    
    void setSources(int n_sources,
                    std::vector<a64::vector<T>> const & deinterlaced_coeffs,
                    SetupParam const & setup,
                    WorkCplxFreqs & work) {
        algo.setup(setup);

        states.clear();
        states.reserve(deinterlaced_coeffs.size());
        for(auto & coeffs : deinterlaced_coeffs) {
            states.push_back(std::make_unique<State>(work));
            auto & c = states.back();
            c->setup(setup);
            c->setCoefficients(std::move(coeffs), algo);
        }
        int const nConvolutionsPerSource = states.size() / n_sources;
        if((n_sources * nConvolutionsPerSource) != states.size()) {
            throw std::runtime_error("inconsistent number of audio sources / channels");
        }
        nConvolutionsPerEar = states.size() / nEars;
        Assert((nEars * nConvolutionsPerEar) == states.size());
        
        int n = 0;
        int const total = states.size();
        for(auto & c : states) {
            dephase(total, n, *c, algo);
            ++n;
        }
    }
    
    void logReport(std::ostream & os) const
    {
        os << "States:" << std::endl;
        IndentingOStreambuf in(os);
        
        int i=0;
        foreachConvReverb([&os, &i, this](auto const & r){
            ++i;
            os << i << ":" << std::endl;
            IndentingOStreambuf indent(os);
            r.logComputeState(os, algo);
        });
    }
    
    int countScales() {
        if (states.empty()) {
            return 0;
        }
        if constexpr(State::has_subsampling) {
            return imajuscule::countScales(*states[0]);
        }
        else {
            return 1;
        }
    }
    
    bool handlesCoefficients() const {
        if(states.empty()) {
            return false;
        }
        return algo.handlesCoefficients();
    }
    
    Latency getLatency() const {
        Assert(handlesCoefficients());
        return states.empty() ? Latency(0) : algo.getLatency();
    }
    
    double getEpsilon() const {
        return epsilonOfNaiveSummation(states, algo) / nEars;
    }
    
    bool isValid() const {
        if(states.empty()) {
            return true;
        }
        return algo.isValid();
    }
    
    template<typename F>
    void foreachConvReverb(F f) const {
        for(auto const & c : states) {
            f(*c);
        }
    }
    
    template<typename F>
    void foreachConvReverb(F f) {
        for(auto & c : states) {
            f(*c);
        }
    }
    
    template<typename FPT2>
    bool assignWetVectorized(FPT2 const * const * const input_buffers,
                             int nInputBuffers,
                             FPT2 ** output_buffers,
                             int nOutputBuffers,
                             int nFramesToCompute,
                             int vectorLength) {
        assert(vectorLength > 0);
        bool success = true;
        
        int i_in = 0;
        auto itConv = states.begin();
        
        Assert(nOutputBuffers == nEars);
        for(int i_out=0; i_out < nOutputBuffers; ++i_out) {
            FPT2 * out = output_buffers[i_out];
            
            bool assign = true;
            for(auto end = i_in+nConvolutionsPerEar;
                i_in < end;
                ++i_in, ++itConv) {
                Assert(i_in < nInputBuffers);
                
                /*std::cout << "out " << i_out << " in " << i_in << std::endl;
                IndentingOStreambuf ind(std::cout);*/

                auto & c = *itConv;
                FPT2 const * const in = input_buffers[i_in];
                for(int i=0; i<nFramesToCompute; i += vectorLength) {
                    if(assign) {
                        c->stepAssignVectorized(in + i,
                                                out + i,
                                                std::min(vectorLength, nFramesToCompute-i),
                                                algo);
                    }
                    else {
                        c->stepAddVectorized(in + i,
                                             out + i,
                                             std::min(vectorLength, nFramesToCompute-i),
                                             algo);
                    }
                }
                assign = false;
                if constexpr (State::step_can_error) {
                    success = !c->hasStepErrors() && success;
                }
            }
            if(unlikely(assign)) {
                fft::RealSignal_<fft::Fastest, FPT2>::zero_n_raw(out, nFramesToCompute);
            }
            if(i_in == nInputBuffers) {
                i_in = 0;
            }
        }
        return success;
    }
    
    template<typename FPT2>
    bool addWetInputZeroVectorized(FPT2 ** output_buffers,
                                   int nOutputBuffers,
                                   int nFramesToCompute,
                                   int vectorLength) {
        assert(vectorLength > 0);
        bool success = true;
        
        Assert(nOutputBuffers == nEars);
        for(int i_out=0; i_out < nOutputBuffers; ++i_out) {
            FPT2 * out = output_buffers[i_out];
            for(int i=0; i<nFramesToCompute; i += vectorLength) {
                for(int i_in = 0; i_in < nConvolutionsPerEar; ++i_in) {
                    auto & c = states[nConvolutionsPerEar*i_out + i_in];
                    c->stepAddInputZeroVectorized(out+i,
                                                  std::min(vectorLength, nFramesToCompute-i),
                                                  algo);
                    if constexpr (State::step_can_error) {
                        success = !c->hasStepErrors() && success;
                    }
                }
            }
        }
        return success;
    }
    
    
    void clear() {
        states.clear();
        nConvolutionsPerEar = 0;
    }
    
    void flushToSilence() {
        for(auto & c : states) {
            c->flushToSilence(algo);
        }
    }
    
private:
    
    int nConvolutionsPerEar = 0;
    std::vector<std::unique_ptr<State>> states;
    Algo algo;

    bool empty() const {
        return nConvolutionsPerEar == 0;
    }
    
};


template<typename SetupParam>
void padForScales(SetupParam const & spec,
                  std::vector<a64::vector<double>> & deinterlaced_coeffs) {
    using namespace std;
    
    int const total_response_size = deinterlaced_coeffs.empty() ? 0 : deinterlaced_coeffs[0].size();
    for(auto const & v:deinterlaced_coeffs) {
        if(total_response_size != v.size()) {
            throw std::runtime_error("deinterlaced coefficients have different sizes");
        }
    }
    
    int total_response_size_padded = 0;
    
    if constexpr (SetupParam::has_subsampling) {
        int const n_scales = count_scales(spec);
        assert(n_scales >= 1);
        int lateHandlerFirstScalePartitionSize = spec.b.a.partition_size;
        int const n_coeffs_early_handler = std::max(minLatencyLateHandlerWhenEarlyHandlerIsDefaultOptimizedFIRFilter,
                                                    earliestDeepestLatency<typename SetupParam::BParam>(lateHandlerFirstScalePartitionSize)).toInteger();
        int const late_response_sz = std::max(0
                                              ,total_response_size - n_coeffs_early_handler);
        int const scale_sz = SameSizeScales::get_scale_sz(late_response_sz, n_scales);
        
        if(n_scales > 1) {
            // pad the coefficients so that all scales have the same rythm.
            int const target_late_response_sz = SameSizeScales::get_max_response_sz(n_scales, scale_sz);
            
            total_response_size_padded = n_coeffs_early_handler + target_late_response_sz;
        }
        else {
            total_response_size_padded = total_response_size;
        }
        
        for(auto & v : deinterlaced_coeffs) {
            v.resize(total_response_size_padded);
        }
    }
    else {
        total_response_size_padded = total_response_size;
    }
}

template<typename Reverb, typename ...Args>
auto findBestParams(int const n_response_channels,
                    int const n_frames,
                    int n_audiocb_frames,
                    double sampleRate,
                    std::ostream & os,
                    Args... args)
{
    using PartitionAlgo = typename Reverb::PartitionAlgo;
    constexpr auto nEars = Reverb::nEars;
    
    if(n_audiocb_frames <= 0) {
        std::stringstream ss;
        ss << "Negative or zero callback size (" << n_audiocb_frames << ")";
        throw std::runtime_error(ss.str());
    }
    
    using namespace std;
    
    static constexpr double ratio_hard_limit = 1.0;
    //    because of overhead due to os managing audio, because of "other things running on the device", etc...
    // at 0.38f on ios release we have glitches when putting the app in the background
    // at 0.25f on linux we have glitches
    static constexpr double ratio_soft_limit = 0.15 * ratio_hard_limit;
    
    double const theoretical_max_ns_per_frame(1e9/sampleRate);
    double const max_avg_time_per_sample(theoretical_max_ns_per_frame * ratio_soft_limit / static_cast<float>(n_response_channels));
    
    auto partitionning = PartitionAlgo::run(n_response_channels,
                                            nEars,
                                            n_audiocb_frames,
                                            n_frames,
                                            sampleRate,
                                            max_avg_time_per_sample,
                                            os,
                                            args...);
    if(!partitionning) {
        std::stringstream ss;
        ss << "could not optimize (2) :" << std::endl << os.rdbuf();
        throw std::runtime_error(ss.str());
    }
    partitionning->logReport(n_response_channels,
                             theoretical_max_ns_per_frame,
                             os);
    return *partitionning;
}

template<typename Reverb>
void applyParams(Reverb & rev,
                 int const n_sources,
                 DeinterlacedBuffers<typename Reverb::FPT> const & deinterlaced,
                 typename Reverb::WorkCplxFreqs & work,
                 typename Reverb::SetupParam const & p,
                 std::ostream & os)
{
    // buffers are copied here:
    auto buffers = deinterlaced.getBuffers();
    padForScales(p, buffers);

    rev.clear();
    rev.setSources(n_sources,
                   buffers,
                   p,
                   work);
    rev.logReport(os);
}

template<typename Reverb, typename ...Args>
void applyBestParams(Reverb & rev,
                     typename Reverb::MemResource::type & memory,
                     int const n_sources,
                     DeinterlacedBuffers<typename Reverb::FPT> const & deinterlaced,
                     typename Reverb::WorkCplxFreqs & work,
                     int n_audiocb_frames,
                     double sampleRate,
                     std::ostream & os,
                     Args... args)
{
    auto p = findBestParams<Reverb>(deinterlaced.countChannels(),
                                    deinterlaced.countFrames(),
                                    n_audiocb_frames,
                                    sampleRate,
                                    os,
                                    args...);
    
    int const sz = Reverb::getAllocationSz(p, deinterlaced.countChannels());
    
    memory.clear();
    memory.reserve(sz);
    
    using MemRsc = typename Reverb::MemResource;
    
    {
        typename MemRsc::use_type use(memory);

        applyParams(rev,
                    n_sources,
                    deinterlaced,
                    work,
                    p,
                    os);
    }
    
    if constexpr (MemRsc::limited) {
        //std::cout << "mem monotonic " << memory.used() << std::endl;
        //std::cout << "mem monotonic " << memory.remaining() << std::endl;
        // verify that getAllocationSz_Setup / getAllocationSz_SetCoefficients are correct
        Assert(sz == memory.used()); // or there is padding to avoid false sharing?
    }
}

static std::ostream & operator << (std::ostream &ss, const ConvReverbOptimizationReportSuccess & r) {
    ss << r.optimizationReport << std::endl;
    return ss;
}
}
