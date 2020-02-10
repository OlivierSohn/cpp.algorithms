
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
    
    using Algo =
    std::conditional_t< reverbType==ReverbType::Offline,
    AlgoOptimizedFIRFilter<FPT, Allocator, Tag>,
    
    std::conditional_t< reverbType==ReverbType::Realtime_Synchronous,
    AlgoZeroLatencyScaledFineGrainedPartitionnedConvolution<FPT, Allocator, Tag>,
    
    //std::conditional_t< reverbType==ReverbType::Realtime_Synchronous_Subsampled,
    //  ZeroLatencyScaledFineGrainedPartitionnedSubsampledConvolution<double, Allocator, Tag>,
    
    std::conditional_t< reverbType==ReverbType::Realtime_Asynchronous_Legacy,
    AlgoZeroLatencyScaledAsyncConvolution<FPT, Allocator, Tag, OnWorkerTooSlow>,
    
    std::conditional_t< reverbType==ReverbType::Realtime_Asynchronous,
    AlgoZeroLatencyScaledAsyncConvolutionOptimized<FPT, Allocator, Tag, OnWorkerTooSlow>,
    
    void
    
    >>>>/*>*/;
    
    using Desc = typename Algo::Desc;
    using State = typename Algo::State;

    using SetupParam = typename Algo::SetupParam;
    using PartitionAlgo = PartitionAlgo<SetupParam, FPT, Tag>;
    using WorkCplxFreqs = typename fft::RealFBins_<Tag, FPT, aP::Alloc>::type;
    
    static int getAllocationSz(SetupParam const & p,
                               int const n_sources,
                               int const n_channels) {
        MinSizeRequirement req = p.getMinSizeRequirement();
        int const input_req = XAndFFTS<FPT, Allocator, Tag>::getAllocationSz_Resize(req.xFftSizes);

        return
        n_sources * input_req +
        n_channels * State::getAllocationSz_SetCoefficients(p);
    }

    // for wir files of wave, it seems the order is by "ears" then by "source":
    // ear Left source 1
    // ...
    // ear Left source nConvolutionsPerEar
    // ear Right source 1
    // ...
    // ear Right source nConvolutionsPerEar
    
    void setSources(int n_sources,
                    std::vector<a64::vector<T>> const & deinterlaced_coeffs,
                    SetupParam const & p,
                    WorkCplxFreqs & work) {
        algo.setup(p);

        input_states.clear();
        
        int const n_channels = deinterlaced_coeffs.size();
        int const nConvolutionsPerSource = n_channels / n_sources;
        if((n_sources * nConvolutionsPerSource) != n_channels) {
            throw std::runtime_error("inconsistent number of audio sources / channels");
        }
        nConvolutionsPerEar = n_channels / nEars;
        Assert((nEars * nConvolutionsPerEar) == n_channels);

        input_states.resize(n_sources);
        
        MinSizeRequirement req = p.getMinSizeRequirement();
        
        work.reserve(req.minWorkSize);
        
        for(auto & i : input_states) {
            auto & x_and_ffts = i.x_and_ffts;
            x_and_ffts.resize(req.minXSize,
                              req.xFftSizes);
            
            x_and_ffts.setPhasePeriod(p);
        }
        
        for(auto & y : outputs) {
            y.resize(req.minYSize);
            // set y progress such that results are written in a single chunk
            if(y.ySz >= 1) {
                y.increment();
            }
        }
        
        auto it = input_states.begin();
        
        for(auto & coeffs : deinterlaced_coeffs) {
            it->channels.push_back(std::make_unique<State>());
            auto & c = it->channels.back();
            c->setCoefficients(algo, std::move(coeffs));

            ++it;
            if(it == input_states.end()) {
                it = input_states.begin();
            }
        }
        
        this->work = &work;

        handlePhases();
    }

    void logReport(std::ostream & os) const
    {
        os << "States:" << std::endl;
        IndentingOStreambuf in(os);
        
        int source_idx=0;
        for(auto const & i : input_states) {
            ++source_idx;
            os << source_idx << ":" << std::endl;
            IndentingOStreambuf in(os);
            
            i.x_and_ffts.logComputeState(os);
            
            int conv_idx = 0;
            for(auto & c : i.channels) {
                ++conv_idx;
                os << conv_idx << ":" << std::endl;
                IndentingOStreambuf indent(os);
                c->logComputeState(algo, os);
            }
        }
    }
    
    int countScales() {
        if (input_states.empty() || input_states[0].channels.empty()) {
            return 0;
        }
        if constexpr(Desc::has_subsampling) {
            return imajuscule::countScales(*(input_states[0].channels[0]));
        }
        else {
            return 1;
        }
    }
    
    bool handlesCoefficients() const {
        if(input_states.empty()) {
            return false;
        }
        return algo.handlesCoefficients();
    }
    
    Latency getLatency() const {
        Assert(handlesCoefficients());
        return input_states.empty() ? Latency(0) : algo.getLatency();
    }
    
    double getEpsilon() const {
        if(input_states.empty() || input_states[0].channels.empty()) {
            return 0;
        }
        return nConvolutionsPerEar * input_states[0].channels[0]->getEpsilon(algo);
    }
    
    bool isValid() const {
        if(input_states.empty()) {
            return true;
        }
        return algo.isValid();
    }
    
    template<typename F>
    void foreachConvReverb(F f) const {
        for(auto const & i : input_states) {
            for(auto const & c:i.channels) {
                f(*c);
            }
        }
    }
    
    template<typename F>
    void foreachConvReverb(F f) {
        for(auto & i : input_states) {
            for(auto & c:i.channels) {
                f(*c);
            }
        }
    }
    
    template<typename FPT2>
    bool assignWetVectorized(FPT2 const * const * const input_buffers,
                             int const nInputBuffers,
                             FPT2 ** output_buffers,
                             int const nOutputBuffers,
                             int const nFramesToCompute,
                             int const vectorLength) {
        assert(vectorLength > 0);
        
        if(unlikely(nConvolutionsPerEar == 0)) {
            for(int o=0; o<nOutputBuffers; ++o) {
                fft::RealSignal_<fft::Fastest, FPT2>::zero_n_raw(output_buffers[o], nFramesToCompute);
            }
            return true;
        }

        Assert(nInputBuffers == input_states.size());
        Assert(nOutputBuffers == outputs.size());

        auto * workData = work->data();

        for(int f=0; f<nFramesToCompute; ++f) {
            int o=0;
            
            for(int i=0; i<nInputBuffers; ++i) {
                auto & input_state = input_states[i];
                input_state.x_and_ffts.push(input_buffers[i][f]);
                
                for(auto & c : input_state.channels) {
                    algo.step(*c,
                              input_state.x_and_ffts,
                              outputs[o],
                              workData);

                    ++o;
                    if(o==nOutputBuffers) {
                        o=0;
                    }
                }
            }
            
            Assert(o == 0);
            
            for(; o<nOutputBuffers; ++o) {
                output_buffers[o][f] = outputs[o].getCurrentSignal();
                outputs[o].increment();
            }
        }
        
        if constexpr (Desc::step_can_error) {
            for(auto const & i : input_states) {
                for(auto const & c : i.channels) {
                    if(c->hasStepErrors()) {
                        return false;
                    }
                }
            }
        }
        
        return true;
    }
    
    template<typename FPT2>
    bool addWetInputZeroVectorized(FPT2 ** output_buffers,
                                   int nOutputBuffers,
                                   int nFramesToCompute,
                                   int vectorLength) {
        assert(vectorLength > 0);
        
        Assert(nOutputBuffers == outputs.size());

        if(unlikely(nConvolutionsPerEar == 0)) {
            return true;
        }
        
        auto * workData = work->data();
        
        // ignore vectorization for now.
        // (To vectorize, we would need to modify x and y so that they have bigger buffers, and modify states)
        for(int f=0; f<nFramesToCompute; ++f) {
            int o=0;
            
            for(auto & input_state : input_states) {
                input_state.x_and_ffts.push(0);
                
                for(auto & c : input_state.channels) {
                    algo.step(*c,
                              input_state.x_and_ffts,
                              outputs[o],
                              workData);

                    ++o;
                    if(o==nOutputBuffers) {
                        o=0;
                    }
                }
            }
            
            Assert(o == 0);
            
            for(; o<nOutputBuffers; ++o) {
                output_buffers[o][f] += outputs[o].getCurrentSignal();
                outputs[o].increment();
            }
        }
        
        if constexpr (Desc::step_can_error) {
            for(auto const & i : input_states) {
                for(auto const & c : i.channels) {
                    if(c->hasStepErrors()) {
                        return false;
                    }
                }
            }
        }
        
        return true;
    }
    
    
    void clear() {
        input_states.clear();
        nConvolutionsPerEar = 0;
    }
    
    void flushToSilence() {
        bool xy_linked = XYLinked();
        int n = 0;

        for(auto & i : input_states) {
            for(auto & state : i.channels) {
                algo.flushToSilence(*state);
            }
            
            int const n_steps = i.x_and_ffts.flushToSilence();
            
            // if a single x corresponds to a single y, we keep the corresponding y in phase with x.
            // else, we keep every y in phase with the first x
            if(xy_linked) {
                outputs[n].flushToSilence(n_steps);
            }
            else if (n==0) {
                for(auto & o : outputs) {
                    o.flushToSilence(n_steps);
                }
            }

            ++n;
        }

        handlePhases();
    }
    
private:
    
    int nConvolutionsPerEar = 0;
    WorkCplxFreqs * work = nullptr;
    Algo algo;
    struct InputAndChannels {
        XAndFFTS<FPT, Allocator, Tag> x_and_ffts;
        std::vector<std::unique_ptr<State>> channels;
    };
    std::vector<InputAndChannels> input_states;
    std::array<Y<FPT, Tag>, nEars> outputs;


    bool empty() const {
        return nConvolutionsPerEar == 0;
    }
    
    bool XYLinked() const {
        if(outputs.size() != input_states.size()) {
            return false;
        }
        for(auto & i : input_states) {
            if(i.channels.size() != 1) {
                return false;
            }
        }
        return true;
    }
    
    void handlePhases() {
        bool xy_linked = XYLinked();
        int n = 0;
        int const total = input_states.size();
        for(auto & i : input_states) {
            float const ratio = n / static_cast<float>(total);
            
            i.x_and_ffts.phase_group_ratio = ratio;
            i.x_and_ffts.dephase([&i, this, xy_linked, n](int x_progress){
                for(auto & state : i.channels) {
                    algo.dephaseStep(*state,
                                     x_progress);
                }

                // if a single x corresponds to a single y, we keep the corresponding y in phase with x.
                // else, we keep every y in phase with the first x
                if(xy_linked) {
                    outputs[n].increment();
                }
                else if (n==0) {
                    for(auto & o : outputs) {
                        o.increment();
                    }
                }
            });

            for(auto & state : i.channels) {
                state->onContextFronteer([ratio](auto & s, auto const & fronteerAlgo){
                    s.dephaseByGroupRatio(ratio, fronteerAlgo);
                });
            }
            ++n;
        }
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
    
    if constexpr (SetupParam::has_subsampling) {
        int const n_scales = count_scales(spec);
        assert(n_scales >= 1);
        int lateHandlerFirstScalePartitionSize = spec.b.a.partition_size;
        int const n_coeffs_early_handler = std::max(minLatencyLateHandlerWhenEarlyHandlerIsDefaultOptimizedFIRFilter,
                                                    earliestDeepestLatency<typename SetupParam::BParam>(lateHandlerFirstScalePartitionSize)).toInteger();
        int const late_response_sz = std::max(0
                                              ,total_response_size - n_coeffs_early_handler);
        int const scale_sz = SameSizeScales::get_scale_sz(late_response_sz, n_scales);
        
        int total_response_size_padded = 0;
        
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
}

template<typename Reverb, typename ...Args>
auto findBestParams(int const n_sources,
                    int const n_response_channels,
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
    
    auto partitionning = PartitionAlgo::run(n_sources,
                                            n_response_channels,
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
    auto p = findBestParams<Reverb>(n_sources,
                                    deinterlaced.countChannels(),
                                    deinterlaced.countFrames(),
                                    n_audiocb_frames,
                                    sampleRate,
                                    os,
                                    args...);
    
    int const sz = Reverb::getAllocationSz(p,
                                           n_sources,
                                           deinterlaced.countChannels());
    
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
