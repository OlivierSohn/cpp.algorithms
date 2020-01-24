
namespace imajuscule {

/*
 
 On filters, and how to chose them (here, filter and convolution means the same thing).
 
 Summary
 -------
 
 ********************************************************************************
 *
 * For ** live ** audio processing, use 'ZeroLatencyScaledFineGrainedPartitionnedConvolution':
 *   it has the smallest worst cost per audio callback.
 *
 * For ** offline ** audio processing, use 'OptimizedFIRFilter':
 *   it has the smallest average cost per sample.
 *
 * Both of them are 0-latency and are designed to scale both for very high and very low count of coefficients.
 ********************************************************************************
 
 Details
 -------
 
 There are 3 metrics that can be optimized:
 
 - average cost :
 How much CPU will I need, on average, to compute a single sample?
 This metric should be minimized when doing offline audio processing, so as to ensure that
 the task is executed as fast as possible.
 
 - "worst audiocallback" cost :
 How much CPU will I need, at worst, during a single audio callback, to compute the corresponding samples?
 Note that the size of the audio callback is important here.
 This metric should be optimized when doing live audio processing,
 to ensure that the audio callback meets its deadline.
 
 - latency :
 does my filter / convolution induce any latency?
 This metric should be optimized when doing live audio processing
 with an instrumentist playing a virtual instrument, because
 it's very hard / unpleasant to play an instrument that has noticeable latency.
 
 Here is the mapping from filter type to the kind of metric that is optimized by that filter:
 
 |----------------------------------------------------|----------------------------|-----------|
 |                                                    |     Optimal time cost      |           |
 |                                                    |----------------------------|           |
 |                                                    | on average | in worst case |           |
 | Filter type                                        | per sample | per callback  | 0-latency |
 |----------------------------------------------------|------------|---------------|-----------|
 | FIRFilter                                          |     .      |       .       |    X      |
 | OptimizedFIRFilter                                 |     X      |       .       |    X      |
 | FinegrainedPartitionnedFFTConvolution              |     .      |       X       |    .      |
 | ZeroLatencyScaledFineGrainedPartitionnedConvolution|     .      |       X       |    X      |
 |----------------------------------------------------|------------|---------------|-----------|
 
 */

template<typename T, template<typename> typename Allocator, typename FFTTag>
using OptimizedFIRFilter =
SplitConvolution <
/**/FIRFilter<T, Allocator>,
/**/CustomScaleConvolution<
/*  */FFTConvolutionIntermediate < PartitionnedFFTConvolutionCRTP<T, Allocator, FFTTag> >>>;

template<typename T, template<typename> typename Allocator, typename FFTTag, LatencySemantic Lat = LatencySemantic::DiracPeak>
using ZeroLatencyScaledFineGrainedPartitionnedSubsampledConvolution =
SplitConvolution<
/**/OptimizedFIRFilter<T, Allocator, FFTTag>,

/**/SplitConvolution<
/*  */FinegrainedPartitionnedFFTConvolution<T, Allocator, FFTTag>, // full resolution
/*  */SubSampled< Lat,
/*    */Delayed<

/**/SplitConvolution<
/*  */FinegrainedPartitionnedFFTConvolution<T, Allocator, FFTTag>, // half resolution
/*  */SubSampled<Lat,
/*    */Delayed<

/**/SplitConvolution<
/*  */FinegrainedPartitionnedFFTConvolution<T, Allocator, FFTTag>, // quarter resolution
/*  */SubSampled<Lat,
/*    */Delayed<

/**/FinegrainedPartitionnedFFTConvolution<T, Allocator, FFTTag> // heighth resolution
>>>>>>>>>>;

template<typename T, template<typename> typename Allocator, typename FFTTag>
using ZeroLatencyScaledFineGrainedPartitionnedConvolution =
SplitConvolution<
/**/OptimizedFIRFilter<T, Allocator, FFTTag>,
/**/FinegrainedPartitionnedFFTConvolution<T, Allocator, FFTTag>>;

template<typename T, template<typename> typename Allocator, typename FFTTag, PolicyOnWorkerTooSlow OnWorkerTooSlow>
using ZeroLatencyScaledAsyncConvolution =
SplitConvolution<
/**/OptimizedFIRFilter<T, Allocator, FFTTag>,
/**/AsyncCPUConvolution <
/*  */CustomScaleConvolution<
/*    */FFTConvolutionIntermediate < PartitionnedFFTConvolutionCRTP<T, Allocator, FFTTag> >>, OnWorkerTooSlow>>;

template<typename T, template<typename> typename Allocator, typename FFTTag, PolicyOnWorkerTooSlow OnWorkerTooSlow>
using ZeroLatencyScaledAsyncConvolutionOptimized =
SplitConvolution<
/**/ZeroLatencyScaledFineGrainedPartitionnedConvolution<T, Allocator, FFTTag>,
/**/AsyncCPUConvolution <
/*  */CustomScaleConvolution<
/*    */FFTConvolutionIntermediate < PartitionnedFFTConvolutionCRTP<T, Allocator, FFTTag> >>, OnWorkerTooSlow>>;


/*
Latency characteristics for each PartitionAlgo:

(non-split leaves)

FIRFilter                             : fixed, 0
FinegrainedPartitionnedFFTConvolution : deduced from optimized partition size

(non-split internals)

CustomScaleConvolution : fixed, deduced from firstSz
AsyncCPUConvolution    : deduced from optimized queueSize

(split)

OptimizedFIRFilter  : latence fixe 0
ZeroLatency***      : latence fixe 0
*/


template<typename A>
struct PartitionAlgo<CustomScaleConvolution<A>> {
    using Convolution = CustomScaleConvolution<A>;
    using SetupParam = typename Convolution::SetupParam;
    using ScalingParam = typename SetupParam::ScalingParam;
    using PS = std::optional<SetupParam>;

    static PS run(int const n_channels,
                  int const n_audio_channels,
                  int const n_audio_frames_per_cb,
                  int const zero_latency_response_size,
                  double const frame_rate,
                  double const max_avg_time_per_sample,
                  std::ostream & os,
                  int const firstSz,
                  XFFtCostFactors const & xFftCostFactors)
    {
        os << "Optimization of CustomScaleConvolution" << std::endl;
        IndentingOStreambuf i(os);

        int const nEarlyCoeffs = std::min(zero_latency_response_size,
                                          firstSz-1);
        int const nLateCoeffs = std::max(0,
                                         zero_latency_response_size - nEarlyCoeffs);

        std::optional<std::pair<std::vector<Scaling>, double>> best =
        getOptimalScalingScheme_ForTotalCpu_ByVirtualSimulation<Convolution>(firstSz,
                                                                             nLateCoeffs,
                                                                             xFftCostFactors);
        PS ps;
        if(best) {
            auto scalingParams = scalingsToParams<ScalingParam>(best->first);
            
            ps = SetupParam{scalingParams};
            ps->setCost(best->second);
        }
        else {
            // dans quels cas arrive-t-on ici?
        }
        return ps;
    }
};

template<typename T, template<typename> typename Allocator, typename FFTTag>
struct PartitionAlgo< OptimizedFIRFilter<T, Allocator, FFTTag> > {
    using Convolution = OptimizedFIRFilter<T, Allocator, FFTTag>;
    using LateHandler = typename Convolution::LateHandler;
    using EarlyHandler = typename Convolution::EarlyHandler;
    using SetupParam = typename Convolution::SetupParam;
    using PS = std::optional<SetupParam>;
    
    static PS run(int n_channels,
                  int n_audio_channels,
                  int n_audio_frames_per_cb,
                  int zero_latency_response_size,
                  double frame_rate,
                  double max_avg_time_per_sample,
                  std::ostream & os,
                  XFFtCostFactors const & xFftCostFactors)
    {
        os << "Optimization of OptimizedFIRFilter" << std::endl;
        IndentingOStreambuf i(os);

        auto const nAsyncScalesDropped = ScaleConvolution_::nDroppedOptimalFor_Split_Bruteforce_Fft;
        int const firstSz = static_cast<int>(pow2(nAsyncScalesDropped.toInteger()));
                
        auto pSpecsLate = PartitionAlgo<LateHandler>::run(n_channels,
                                                          n_audio_channels,
                                                          n_audio_frames_per_cb,
                                                          zero_latency_response_size,
                                                          frame_rate,
                                                          max_avg_time_per_sample,
                                                          os,
                                                          firstSz,
                                                          xFftCostFactors
                                                          );
        PS ps;
        if(pSpecsLate) {
            ps = {
                {
                    pSpecsLate->handlesCoefficients() ?
                    pSpecsLate->getImpliedLatency().toInteger() :
                    zero_latency_response_size
                },
                *pSpecsLate
            };
            ps->setCost(pSpecsLate->getCost());
        }
        return ps;
    }
};

enum class SimulationPhasingMode {
    On,
    Off
};
struct SimulationPhasing {
    SimulationPhasingMode mode;
    std::optional<int> groupSize;
    
    static SimulationPhasing no_phasing() {
        return {SimulationPhasingMode::Off, {}};
    }
    static SimulationPhasing phasing_with_group_size(int sz) {
        return {SimulationPhasingMode::On, {sz}};
    }
};


template<typename F>
int computeQueueSize(F nextProcessingDuration,
                     double const producerPeriod,
                     int nIterations) {
    double time = 0;
    int consecutiveMiss = 0;
    int minQueueSz = 1;
    double timeNextProduced = time + producerPeriod;
    double timeNextProcessed = time + nextProcessingDuration();
    while(nIterations >= 0) {
        if(timeNextProduced < timeNextProcessed) {
            ++consecutiveMiss;
            timeNextProduced += producerPeriod;
        }
        else {
            if(consecutiveMiss) {
                --consecutiveMiss;
            }
            timeNextProcessed += nextProcessingDuration();
            --nIterations;
        }
        Assert(consecutiveMiss >= 0);
        minQueueSz = std::max(minQueueSz, consecutiveMiss);
    }
    return minQueueSz;
}


template<typename Algo, PolicyOnWorkerTooSlow OnWorkerTooSlow>
struct PartitionAlgo< AsyncCPUConvolution<Algo, OnWorkerTooSlow> > {
    using Convolution = AsyncCPUConvolution<Algo, OnWorkerTooSlow>;
    using AsyncPart = typename Convolution::AsyncPart;
    using AsyncSetupParam = typename AsyncPart::SetupParam;
    using SetupParam = typename Convolution::SetupParam;
    using PS = std::optional<SetupParam>;
    
    // to account for a system that is under heavier load
    static constexpr int safeFactor = 4;
    
    static PS run(int const n_channels,
                  int const n_audio_channels,
                  int const n_audio_frames_per_cb,
                  int const zero_latency_response_size,
                  double const frame_rate,
                  double const max_avg_time_per_sample,
                  std::ostream & os,
                  SimulationPhasing const & phasing,
                  CountDroppedScales const & nAsyncScalesDropped)
    {
        os << "Optimization of AsyncCPUConvolution" << std::endl;
        IndentingOStreambuf i(os);

        if(n_channels <= 0) {
            throw std::runtime_error("0 channel");
        }
        if(n_audio_frames_per_cb <= 0) {
            throw std::runtime_error("0 n_audio_frames_per_cb");
        }
        int const submission_period = n_audio_frames_per_cb;
        int const nMinDroppedCoeffs = submission_period; // because final latency will be at least submission_period
        int const nAsyncCoeffs = zero_latency_response_size-nMinDroppedCoeffs;
        std::optional<int> maybeQueueSz;
        std::optional<std::pair<std::vector<Scaling>, double>> optimalScaling;
                
        AsyncSetupParam scalingParams({});

        if(nAsyncCoeffs <= 0) {
            maybeQueueSz = 0;
        }
        else {
            int const firstSz = static_cast<int>(pow2(nAsyncScalesDropped.toInteger()));
            XFFtCostFactors unbiasedXFftCostFactors;

            // Because of the queueing mechanism,
            // the real number of late coeffs will be _less_ than what this is optimal for.
            // Hence, as soon as we know the size of the queue, we will be able to deduce the exact
            // number of late coefficients and adapt the scaling scheme accordingly.
            //
            // This future adaptation will be minor (we may remove work but _not_ change the structure
            // because we would have to run the simulation again, as the peak computation times could be worse.
            auto scalingParams2 = PartitionAlgo<AsyncPart>::run(n_channels,
                                                                n_audio_channels,
                                                                n_audio_frames_per_cb,
                                                                // so that the actual number of ceofficients handled is nAsyncCoeffs:
                                                                nAsyncCoeffs + firstSz-1,
                                                                frame_rate,
                                                                max_avg_time_per_sample,
                                                                os,
                                                                firstSz,
                                                                unbiasedXFftCostFactors);
            if(!scalingParams2) {
                throw std::runtime_error("no best");
            }
            
            scalingParams = *scalingParams2;
            
            maybeQueueSz = computeQueueMetricsWithVirtualSimulation(n_channels,
                                                                    n_audio_channels,
                                                                    n_audio_frames_per_cb,
                                                                    nAsyncCoeffs,
                                                                    frame_rate,
                                                                    scalingParams,
                                                                    submission_period,
                                                                    phasing,
                                                                    os);
            if(!maybeQueueSz) {
                os << "Synchronous thread is too slow" << std::endl;
                return {};
            }
        }
        
        int const queueSize = *maybeQueueSz;
        os << "queueSize=" << queueSize << std::endl;
        
        SetupParam ps(
            submission_period,
            queueSize,
            scalingParams
        );
        ps.setCost(scalingParams.getCost());
        
        // now that we know what latency we have, we can adjust:

        if(ps.handlesCoefficients()) {
            int const adjustedNAsyncCoeffs = zero_latency_response_size - ps.getImpliedLatency().toInteger();
            Assert(adjustedNAsyncCoeffs <= nAsyncCoeffs);
            ps.adjustWork(adjustedNAsyncCoeffs);
        }
        return ps;
    }
private:
    
    static std::optional<int>
    computeQueueMetricsWithVirtualSimulation(int const n_channels,
                                             int const n_audio_channels,
                                             int const n_audio_frames_per_cb,
                                             int const nLateCoeffs,
                                             double const frame_rate,
                                             AsyncSetupParam const & asyncParam,
                                             int const submission_period_frames,
                                             SimulationPhasing const & phasing,
                                             std::ostream & os)
    {
        // at the moment, every AsyncCPUConvolution has its own worker thread, so we simulate
        // 'n_channels' threads runing at the same time, to see what percentage of the time they
        // are actually running on the cpu, and what their max 'sleep' time is.
        
        auto const nThreads = n_channels;
        
        os << "Virtual simulation with " << nThreads << " threads:" << std::endl;
        IndentingOStreambuf i(os);

        using Sim = typename AsyncPart::Simulation;
        
        std::vector<std::unique_ptr<Sim>> simu_convs;
        simu_convs.reserve(n_channels);
        for(int i=0; i<n_channels; ++i) {
            simu_convs.push_back(std::make_unique<Sim>());
        }
        
        for(auto & c : simu_convs) {
            c->setup(asyncParam);
            c->setCoefficientsCount(nLateCoeffs);
        }
        if(phasing.mode == SimulationPhasingMode::On)
        {
            int n=0;
            int const total_sz = phasing.groupSize ? *phasing.groupSize : simu_convs.size();
            Assert(total_sz);
            for(auto & c : simu_convs) {
                dephase(total_sz,
                        n % total_sz,
                        *c);
                ++n;
            }
        }
                
        MaxWallTimeIncrementEval::BythreadMaxIncrements maxIncrements;
        maxIncrements.resize(nThreads);
        
        // We simulate without taking into account the effects of the realtime audio thread.
        // We will take this into account later.
        
        auto allocs = computeAllocationFactors(nThreads,
                                               std::chrono::seconds(1),
                                               MaxWallTimeIncrementEval(maxIncrements));
        
        std::vector<std::chrono::steady_clock::duration> maxIncrementsValues;
        maxIncrementsValues.resize(maxIncrements.size());
        std::transform(maxIncrements.begin(),
                       maxIncrements.end(),
                       maxIncrementsValues.begin(),
                       [](auto & o) {
            if(!o.value) {
                throw std::runtime_error("no max increment");
            }
            return *o.value;
        });
        
        auto const max_sleep_time = *std::max_element(maxIncrementsValues.begin(),
                                                      maxIncrementsValues.end());
        // we apply the safe factor
        double max_sleep_time_seconds = safeFactor * std::chrono::duration_cast<std::chrono::microseconds>(max_sleep_time).count() / 1000000.;
        
        double const cb_period_seconds = n_audio_frames_per_cb / frame_rate;
        max_sleep_time_seconds = std::max(max_sleep_time_seconds,
                                          // 1 because the audio thread could preempt and run at most cb_period
                                          // 1 to avoid audio dropouts in case of the occurence of the race condition
                                          //   commented in the code of the worker of AsyncCPUConvolution
                                          2. * cb_period_seconds
                                          );
        
        os << "max_sleep_time_seconds " << max_sleep_time_seconds << std::endl;
        
        int const max_sleep_frames = static_cast<int>(0.5 + std::ceil(max_sleep_time_seconds * frame_rate));
        int const max_sleep_submissions = 1 + (max_sleep_frames-1) / n_audio_frames_per_cb;
        Assert(max_sleep_submissions >= 0);
        
        double min_allocation_factor = *std::min_element(allocs.begin(), allocs.end());
        
        // the realtime audio thread also uses some cpu, here we assume the worst case could be that
        // it takes half the cpu - which is a lot.
        min_allocation_factor /= 2.;
        // and we apply the safe factor
        min_allocation_factor /= safeFactor;
        Assert(min_allocation_factor > 0.);
        Assert(min_allocation_factor < 1.);
        
        os << "min_allocation_factor " << min_allocation_factor << std::endl;
        
        double const submission_period_seconds = submission_period_frames / frame_rate;
        
        XFFtCostFactors unbiasedXFftCostFactors;
        int queueSize = computeQueueSize([&simu_convs,
                                          min_allocation_factor,
                                          &unbiasedXFftCostFactors](){
            // take only the first one into account (each of them is equivalent over a period)
            return simu_convs[0]->simuStep(unbiasedXFftCostFactors) / min_allocation_factor;
        },
                                         submission_period_seconds,
                                         2*nLateCoeffs);
        
        os << "queueSize from allocations: " << queueSize << " + from sleep: " << max_sleep_submissions<< std::endl;
        
        return queueSize + max_sleep_submissions;
    }
    
    static std::optional<int>
    computeQueueMetricsWithRealSimulation(int n_channels,
                                          int n_audio_channels,
                                          int n_audio_frames_per_cb,
                                          int nLateCoeffs,
                                          double frame_rate,
                                          std::vector<Scaling> const & scaling,
                                          int const submission_period_frames,
                                          SimulationPhasing const & phasing,
                                          std::ostream & os)
    {
        struct Metrics {
            range<int> resultQueueEltsCountRange;
            
            // We are interested in the minimum size, by period, of the result queue.
            // If the minimum size keeps being smaller over the periods, it meens that
            // the background thread is too slow.
            
            double minResultQueueEltsCountRange_gradient = 0; // diff between last and first period over number of periods
            // a strictly negative gradient means that the system is not fast enough.
            
            int minResultQueueEltsCountRange_min_local_gradient = 0; // diff between consecutive periods
            
            void mergeWorst(Metrics const & o) {
                resultQueueEltsCountRange.extend(o.resultQueueEltsCountRange);
                minResultQueueEltsCountRange_gradient = std::min(minResultQueueEltsCountRange_gradient,
                                                                 o.minResultQueueEltsCountRange_gradient);
                minResultQueueEltsCountRange_min_local_gradient = std::min(minResultQueueEltsCountRange_min_local_gradient,
                                                                           o.minResultQueueEltsCountRange_min_local_gradient);
            }
        };
        
        struct QueueMetrics {
            QueueMetrics(int period, int nPeriods)
            : nPeriods(nPeriods)
            , period_frames(period)
            {
                if(nPeriods < 2) {
                    throw std::logic_error("2 periods are required to compute a gradient");
                }
                resultQueueEltsCountRangeByPeriod.reserve(nPeriods);
                resultQueueEltsCountRangeByPeriod.resize(1);
                signalQueueEltsCountRangeByPeriod.reserve(nPeriods);
                signalQueueEltsCountRangeByPeriod.resize(1);
            }
            
        private:
            int const nPeriods;
            int const period_frames;
            int period = 0;
            int n_cur_frame = 0;
            bool hasErrors = false;
            
            std::vector<range<int>> resultQueueEltsCountRangeByPeriod, signalQueueEltsCountRangeByPeriod;
        public:
            bool recordQueueSize(int const resultQueueEltsCount,
                                 int const signalQueueEltsCount,
                                 bool const hasErrorWorkerTooSlow,
                                 int const nFrames) {
                assert(period < nPeriods);
                if(hasErrorWorkerTooSlow) {
                    hasErrors = true;
                    return false;
                }
                n_cur_frame += nFrames;
                if(n_cur_frame >= period_frames) {
                    n_cur_frame -= period_frames;
                    ++period;
                    if(period >= nPeriods) {
                        return false;
                    }
                    Assert(resultQueueEltsCountRangeByPeriod.size() < resultQueueEltsCountRangeByPeriod.capacity());
                    Assert(signalQueueEltsCountRangeByPeriod.size() < signalQueueEltsCountRangeByPeriod.capacity());
                    resultQueueEltsCountRangeByPeriod.emplace_back();
                    signalQueueEltsCountRangeByPeriod.emplace_back();
                }
                resultQueueEltsCountRangeByPeriod.back().extend(resultQueueEltsCount);
                signalQueueEltsCountRangeByPeriod.back().extend(signalQueueEltsCount);
                return true;
            }
            
            std::vector<range<int>> const & getResultQueueEltsCountRangeByPeriod () const {
                return resultQueueEltsCountRangeByPeriod;
            }
            std::vector<range<int>> const & getSignalQueueEltsCountRangeByPeriod () const {
                return signalQueueEltsCountRangeByPeriod;
            }
            
            std::optional<Metrics> getMetrics() const {
                if(hasErrors) {
                    return {};
                }
                Metrics m;
                for(auto const & r : resultQueueEltsCountRangeByPeriod) {
                    m.resultQueueEltsCountRange.extend(r);
                }
                {
                    std::optional<int> prevMinQueueSize;
                    std::optional<int> minGradient;
                    bool first = true;
                    for(auto const & r : resultQueueEltsCountRangeByPeriod) {
                        if(first) {
                            // discard first period for gradient computation
                            // (we are not in a stable mode at the beginning of the first period)
                            first = false;
                            break;
                        }
                        auto minResultQueueSizeDuringPeriod = r.getMin();
                        if(!prevMinQueueSize) {
                            prevMinQueueSize = minResultQueueSizeDuringPeriod;
                        }
                        else {
                            auto gradient = minResultQueueSizeDuringPeriod - *prevMinQueueSize;
                            if(!minGradient || *minGradient > gradient) {
                                minGradient = gradient;
                            }
                        }
                    }
                    m.minResultQueueEltsCountRange_min_local_gradient = *minGradient;
                }
                if(resultQueueEltsCountRangeByPeriod.size() >= 3)
                {
                    auto & first = *(resultQueueEltsCountRangeByPeriod.begin() + 1);
                    auto & last = resultQueueEltsCountRangeByPeriod.back();
                    m.minResultQueueEltsCountRange_gradient = (last.getMin() - first.getMin()) / static_cast<double>(nPeriods);
                }
                else {
                    throw std::logic_error("not enough periods to compute a gradient");
                }
                return m;
            }
        };
        
        std::vector<std::unique_ptr<Convolution>> async_convs;
        async_convs.reserve(n_channels);
        for(int i=0; i<n_channels; ++i) {
            async_convs.push_back(std::make_unique<Convolution>());
        }
        
        a64::vector<double> coeffs;
        coeffs.resize(nLateCoeffs);
        
        // we use a huge queue so that we can have room for queue size variations
        int const num_frames_in_queue = std::max(10000, nLateCoeffs/2);
        
        int queue_size = num_frames_in_queue / submission_period_frames;
        if(queue_size * submission_period_frames < num_frames_in_queue) {
            ++queue_size;
        }
        
        using ScalingParam = typename AsyncPart::SetupParam::ScalingParam;
        auto scalingParams = scalingsToParams<ScalingParam>(scaling);
        
        for(auto & c : async_convs) {
            c->setup(
                     {
                submission_period_frames,
                queue_size,
                {
                    scalingParams
                }
            });
            c->setCoefficients(coeffs);
        }
        if(phasing.mode == SimulationPhasingMode::On)
        {
            int n=0;
            int const total_sz = phasing.groupSize ? *phasing.groupSize : async_convs.size();
            Assert(total_sz);
            for(auto & c : async_convs) {
                dephase(total_sz,
                        n % total_sz,
                        *c);
                ++n;
            }
        }
        
        AudioHostSimulator simu{
            frame_rate,
            n_audio_frames_per_cb,
            n_channels,
            n_audio_channels
        };
        
        std::vector<QueueMetrics> m;
        m.reserve(async_convs.size());
        for(auto const & c : async_convs) {
            int const metric_period = c->getAsyncAlgo().getBiggestScale();
            
            static constexpr int numPeriods = 5; // the first period is not taken into account for gradient computation
            
            m.emplace_back(metric_period,
                           numPeriods);
        }
        
        auto p = simu.simulate([&async_convs, &m]
                               (std::vector<a64::vector<float>> const & inputs,
                                std::vector<a64::vector<float>> & outputs,
                                int const cur_frame,
                                int const n_frames) {
            Assert(inputs.size() == async_convs.size());
            Assert(!outputs.empty());
            
            for(int i=0; i < async_convs.size(); ++i) {
                auto & conv = *async_convs[i];
                auto const & input = inputs[i];
                auto & output = outputs[i%outputs.size()];
                
                bool const assign = i < outputs.size();
                
                if(assign) {
                    conv.stepAssignVectorized(input.data()+cur_frame,
                                              output.data()+cur_frame,
                                              n_frames);
                }
                else {
                    conv.stepAddVectorized(input.data()+cur_frame,
                                           output.data()+cur_frame,
                                           n_frames);
                }
                
                
                if(!m[i].recordQueueSize(conv.getResultQueueSize(),
                                         conv.getSignalQueueSize(),
                                         conv.hasStepErrors(),
                                         n_frames)) {
                    // the number of periods has elapsed or the worker is too slow
                    return false;
                }
            }
            return true;
        });
        
        if(!p) {
            // missed deadline
            return {};
        }
        auto & periodically = *p;
        auto & queueMetrics = m;
        if(queueMetrics.empty()) {
            // no late coeffs
            return 0;
        }
        else {
            std::vector<Metrics> metrics;
            metrics.reserve(queueMetrics.size());
            for(auto const & m : queueMetrics) {
                auto mayMetric = m.getMetrics();
                if(!mayMetric) {
                    os << "Background worker is too slow : a result queue was empty" << std::endl;
                    return {};
                }
                metrics.push_back(*mayMetric);
            }
            
            Metrics worst;
            for(auto const & m : metrics) {
                worst.mergeWorst(m);
            }
            
            auto normalizedGradient = worst.minResultQueueEltsCountRange_gradient / (1 + worst.resultQueueEltsCountRange.getSpan());
            
            // todo use a better approach than this threshold
            static constexpr double normalizedGradientThreshold = -0.2;
            
            if(normalizedGradient < normalizedGradientThreshold) {
                os << "Background worker is too slow : normalizedGradient on min queue elts count = " << normalizedGradient << std::endl;
                int i=0;
                for(auto const & m : queueMetrics) {
                    std::cout << i << std::endl;
                    for(auto const & r : m.getResultQueueEltsCountRangeByPeriod()) {
                        std::cout << r.getMin() << " " << r.getMax() << std::endl;
                    }
                    std::cout << std::endl;
                    for(auto const & r : m.getSignalQueueEltsCountRangeByPeriod()) {
                        std::cout << r.getMin() << " " << r.getMax() << std::endl;
                    }
                    std::cout << std::endl;
                    ++i;
                }
                return {};
            }
            
            int minQueueSize = // assumes that the submission_period was n_audio_frames_per_cb during the simulation
            1 +
            worst.resultQueueEltsCountRange.getSpan();
            
            // to avoid audio dropouts in case of the occurence of the race condition
            // commented in the code of the worker of AsyncCPUConvolution:
            ++minQueueSize;
            
            minQueueSize *= safeFactor;
            return minQueueSize;
        }
    }
};

template<typename T, template<typename> typename Allocator, typename FFTTag, PolicyOnWorkerTooSlow OnWorkerTooSlow>
struct PartitionAlgo< ZeroLatencyScaledAsyncConvolution<T, Allocator, FFTTag, OnWorkerTooSlow> > {
    using Convolution = ZeroLatencyScaledAsyncConvolution<T, Allocator, FFTTag, OnWorkerTooSlow>;
    using SetupParam = typename Convolution::SetupParam;
    using PS = std::optional<SetupParam>;
    
    using LateHandler = typename Convolution::LateHandler;
    using EarlyHandler = typename Convolution::EarlyHandler;
    
    static PS run(int n_channels,
                  int n_audio_channels,
                  int n_audio_frames_per_cb,
                  int zero_latency_response_size,
                  double frame_rate,
                  double max_avg_time_per_sample,
                  std::ostream & os,
                  SimulationPhasing const & phasing)
    {
        os << "Optimization of ZeroLatencyScaledAsyncConvolution" << std::endl;
        IndentingOStreambuf i(os);

        // There is a balance to find between the risk of audio dropouts due to:
        //   - A : the earlyhandler having too many coefficients to handle
        //   - B : the latehandler having too many (small) scales to handle
        //
        // in the choice of nAsyncScalesDropped, here we chose to care about 'B':
        //
        // We set nAsyncScalesDropped such that, in the async part, we can use scale convolution
        // on the range of ffts where it's more optimal than brute force convolution.
        //
        // But using nDroppedOptimalFor_Split_Bruteforce_Fft instead of 0 for nAsyncScalesDropped
        // incurs more latency in the latehandler, hence more coefficients to handle in the early handler.
        //
        // So this design choice may be changed in the future.
        auto const nAsyncScalesDropped = ScaleConvolution_::nDroppedOptimalFor_Split_Bruteforce_Fft;
        
        auto pSpecsLate = PartitionAlgo<LateHandler>::run(n_channels,
                                                          n_audio_channels,
                                                          n_audio_frames_per_cb,
                                                          zero_latency_response_size,
                                                          frame_rate,
                                                          max_avg_time_per_sample,
                                                          os,
                                                          phasing,
                                                          nAsyncScalesDropped);
        PS ps;
        if(pSpecsLate) {
            int nEarlyCoeffs = pSpecsLate->handlesCoefficients() ? pSpecsLate->getImpliedLatency().toInteger() : zero_latency_response_size;
            XFFtCostFactors unbiasedXFftCostFactors;
            auto pSpecsEarly = PartitionAlgo<EarlyHandler>::run(n_channels,
                                                                n_audio_channels,
                                                                n_audio_frames_per_cb,
                                                                nEarlyCoeffs,
                                                                frame_rate,
                                                                max_avg_time_per_sample,
                                                                os,
                                                                unbiasedXFftCostFactors);
            if(pSpecsEarly) {
                ps = {
                    *pSpecsEarly,
                    *pSpecsLate
                };
                ps->setCost(pSpecsEarly->getCost() +
                            pSpecsLate->getCost());
            }
        }
        return ps;
    }
};

template<typename T, template<typename> typename Allocator, typename FFTTag, PolicyOnWorkerTooSlow OnWorkerTooSlow>
struct PartitionAlgo< ZeroLatencyScaledAsyncConvolutionOptimized<T, Allocator, FFTTag, OnWorkerTooSlow> > {
    using Convolution = ZeroLatencyScaledAsyncConvolutionOptimized<T, Allocator, FFTTag, OnWorkerTooSlow>;
    using SetupParam = typename Convolution::SetupParam;
    using PS = std::optional<SetupParam>;
    
    using EarlyHandler = typename Convolution::EarlyHandler;
    using LateHandler = typename Convolution::LateHandler;
    
    using LateSetupParam = typename LateHandler::SetupParam;
    
    static PS run(int n_channels,
                  int n_audio_channels,
                  int n_audio_frames_per_cb,
                  int zero_latency_response_size,
                  double frame_rate,
                  double max_avg_time_per_sample,
                  std::ostream & os,
                  SimulationPhasing const & phasing)
    {
        os << "Optimization of ZeroLatencyScaledAsyncConvolutionOptimized" << std::endl;
        IndentingOStreambuf i(os);

        // There is a balance to find between the risk of audio dropouts due to:
        //   - A : the earlyhandler having too many coefficients to handle
        //   - B : the latehandler having too many (small) scales to handle
        //
        // in the choice of nAsyncScalesDropped, here we chose to care about 'B':
        //
        // We set nAsyncScalesDropped such that, in the async part, we can use scale convolution
        // on the range of ffts where it's more optimal than brute force convolution.
        //
        // But using nDroppedOptimalFor_Split_Bruteforce_Fft instead of 0 for nAsyncScalesDropped
        // incurs more latency in the latehandler, hence more coefficients to handle in the early handler.
        //
        // So this design choice may be changed in the future.
        auto const nAsyncScalesDropped = ScaleConvolution_::nDroppedOptimalFor_Split_Bruteforce_Fft;
        
        auto pSpecsLate = PartitionAlgo<LateHandler>::run(n_channels,
                                                          n_audio_channels,
                                                          n_audio_frames_per_cb,
                                                          zero_latency_response_size,
                                                          frame_rate,
                                                          max_avg_time_per_sample,
                                                          os,
                                                          phasing,
                                                          nAsyncScalesDropped);
        PS ps;
        if(pSpecsLate) {
            int nEarlyCoeffs = pSpecsLate->handlesCoefficients() ? pSpecsLate->getImpliedLatency().toInteger() : zero_latency_response_size;
            os << "Deduced early handler coeffs:" << nEarlyCoeffs << std::endl;
            
            auto pSpecsEarly = PartitionAlgo<EarlyHandler>::run(n_channels,
                                                                n_audio_channels,
                                                                n_audio_frames_per_cb,
                                                                nEarlyCoeffs,
                                                                frame_rate,
                                                                max_avg_time_per_sample,
                                                                os);
            
            if(pSpecsEarly) {
                ps = {
                    *pSpecsEarly,
                    *pSpecsLate
                };
                ps->setCost(pSpecsEarly->getCost() +
                            pSpecsLate->getCost());
            }
        }
        return ps;
    }
};

  

template<typename T>
struct EarlyestDeepest {
    using type = T;
};

template<typename A, typename B>
struct EarlyestDeepest< SplitConvolution<A,B> > {
    using type = typename EarlyestDeepest<A>::type;
};


template<typename C>
Latency constexpr earliestDeepestLatency(int partition_sz) {
    using ED = typename EarlyestDeepest<C>::type;
    return ED::getLatencyForPartitionSize(partition_sz);
}


// Subsampling can be used to diminish the resolution of the impulse response tail,
// it makes the computations use less CPU cycles:
enum class ResponseTailSubsampling {
    // response is used at full resolution everywhere (most CPU intensive):
    FullRes, // 0
    ScaleCount_1 = FullRes,
    // the beginning of the response is at full resolution, then half resolution:
    UpToHalfRes, // 1
    ScaleCount_2 = UpToHalfRes,
    // the beginning of the response is at full resolution, then half resolution, then quarter resolution:
    UpToQuarterRes, // 2
    ScaleCount_3 = UpToQuarterRes,
    // the beginning of the response is at full resolution, then half resolution, then quarter resolution, then heighth resolution:
    UpToHeighthRes, // 3
    ScaleCount_4 = UpToHeighthRes,
    // If in doubt, use this mode: the least number of scales will be used
    // if we can afford the induced computations (using an auto optimizing algorithm).
    HighestAffordableResolution, // 4
};

template<typename Convolution>
range<int> getScaleCountRanges(ResponseTailSubsampling rts) {
    range<int> r;
    
    if constexpr (Convolution::has_subsampling) {
        switch(rts) {
            case ResponseTailSubsampling::ScaleCount_1:
                r.extend(1);
                break;
            case ResponseTailSubsampling::ScaleCount_2:
                r.extend(2);
                break;
            case ResponseTailSubsampling::ScaleCount_3:
                r.extend(3);
                break;
            case ResponseTailSubsampling::ScaleCount_4:
                r.extend(4);
                break;
            default:
                assert(0);
            case ResponseTailSubsampling::HighestAffordableResolution:
                r.set(1, nMaxScales);
                break;
        }
    }
    else {
        r.extend(1);
    }
    
    return r;
}


template<typename T, template<typename> typename Allocator, typename FFTTag, typename SPEarly, typename SPLate>
auto mkSubsamplingSetupParams(SPEarly const & early_params,
                              SPLate const & late_params,
                              int const n_scales,
                              int const scale_sz) {
    using SetupParam = typename ZeroLatencyScaledFineGrainedPartitionnedSubsampledConvolution<T, Allocator, FFTTag>::SetupParam;
    int delay = 0;
    if(n_scales > 1) {
        delay = SameSizeScales::getDelays(scale_sz,
                                          late_params.getImpliedLatency());
        Assert(delay > 0);
    }
    
    // set the delays
    
    // we disable the unused scales by setting the partition size to 0.
    auto zero = SPLate::makeInactive();
    static_assert(4==nMaxScales);
    std::array<SPLate, nMaxScales> ps {
        late_params,
        (n_scales >= 2)?late_params:zero,
        (n_scales >= 3)?late_params:zero,
        (n_scales >= 4)?late_params:zero
    };
    Assert(n_scales <= nMaxScales);
    return SetupParam
    {
        early_params,
        {
            ps[0],
            {
                {
                    (n_scales >= 2)?delay:0,
                    {
                        ps[1],
                        {
                            {
                                (n_scales >= 3)?delay:0,
                                {
                                    ps[2],
                                    {
                                        {
                                            (n_scales >= 4)?delay:0,
                                            ps[3]
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    };
}

// TODO make more generic
template<typename SetupParam>
int count_scales(SetupParam const & p) {
    int n_scales = 1;
    if(p.b.a.handlesCoefficients()) {
        if(p.b.b.subsampled.delayed.a.handlesCoefficients()) {
            n_scales = 2;
            if(p.b.b.subsampled.delayed.b.subsampled.delayed.a.handlesCoefficients()) {
                n_scales = 3;
                if(p.b.b.subsampled.delayed.b.subsampled.delayed.b.subsampled.delayed.handlesCoefficients()) {
                    n_scales = 4;
                }
            }
        }
    }
    return n_scales;
}

template<typename T, template<typename> typename Allocator, typename FFTTag>
struct PartitionAlgo< ZeroLatencyScaledFineGrainedPartitionnedConvolution<T, Allocator, FFTTag> > {
    using Convolution = ZeroLatencyScaledFineGrainedPartitionnedConvolution<T, Allocator, FFTTag>;
    
private:
    using LateHandler = typename Convolution::LateHandler;
    using EarlyHandler = typename Convolution::EarlyHandler;
public:
    using SetupParam = typename Convolution::SetupParam;
    using PS = std::optional<SetupParam>;
    
    static PS run(int n_channels,
                  int n_audio_channels,
                  int n_audio_frames_per_cb,
                  int zero_latency_response_size,
                  double frame_rate,
                  double max_avg_time_per_sample,
                  std::ostream & os)
    {
        os << "Optimization of ZeroLatencyScaledFineGrainedPartitionnedConvolution" << std::endl;
        IndentingOStreambuf i(os);

        auto getEarlyHandlerParams = [&](int countEarlyHandlerCoeffs,
                                         XFFtCostFactors const & xFftCostFactors) {
            return PartitionAlgo<EarlyHandler>::run(n_channels,
                                                    n_audio_channels,
                                                    n_audio_frames_per_cb,
                                                    countEarlyHandlerCoeffs,
                                                    frame_rate,
                                                    max_avg_time_per_sample,
                                                    os,
                                                    xFftCostFactors);
            
        };
        if(zero_latency_response_size <= minLatencyLateHandlerWhenEarlyHandlerIsDefaultOptimizedFIRFilter.toInteger()) {
            // in that case there is no late handling at all.
            PS ps;
            XFFtCostFactors unbiasedXFftCostFactors;
            auto earlyRes = getEarlyHandlerParams(zero_latency_response_size,
                                                  unbiasedXFftCostFactors);
            if(earlyRes)
            {
                ps = {
                    *earlyRes,
                    FinegrainedSetupParam::makeInactive()
                };
                ps->setCost(earlyRes->getCost());
            }
            return ps;
        }
        auto lateRes = PartitionAlgo<LateHandler>::run(n_channels,
                                                       1,
                                                       n_audio_frames_per_cb,
                                                       zero_latency_response_size,
                                                       os);
        PS ps;
        if(lateRes) {
            // this fft comes "for free" because the latehandler uses it:
            XFFtCostFactors xFftCostFactors{{
                {lateRes->partition_size, 0.f}
            }};
            auto earlyRes = getEarlyHandlerParams(LateHandler::nCoefficientsFadeIn +
                                                  LateHandler::getLatencyForPartitionSize(lateRes->partition_size).toInteger(),
                                                  xFftCostFactors);
            if(earlyRes) {
                ps = SetupParam {
                    *earlyRes,
                    *lateRes
                };
                ps->setCost(earlyRes->getCost() +
                            lateRes->getCost());
            }
        }
        return ps;
    }
};

template<typename T, template<typename> typename Allocator, typename FFTTag>
struct PartitionAlgo< ZeroLatencyScaledFineGrainedPartitionnedSubsampledConvolution<T, Allocator, FFTTag> > {
    using Convolution = ZeroLatencyScaledFineGrainedPartitionnedSubsampledConvolution<T, Allocator, FFTTag>;
    
private:
    using EarlyHandler = typename Convolution::EarlyHandler;
    using EarliestDeepesLateHandler = typename EarlyestDeepest<typename Convolution::LateHandler>::type;
public:
    using SetupParam = typename Convolution::SetupParam;
    using PS = std::optional<SetupParam>;
    
    static PS run(int const n_response_channels,
                  int const n_audio_channels,
                  int const n_audio_frames_per_cb,
                  int const zero_latency_response_size,
                  double const frame_rate,
                  double const max_avg_time_per_sample,
                  std::ostream & os,
                  ResponseTailSubsampling rts)
    {
        os << "Optimization of ZeroLatencyScaledFineGrainedPartitionnedSubsampledConvolution" << std::endl;
        IndentingOStreambuf i(os);

        PS res;
        
        range<int> const scales = getScaleCountRanges<Convolution>(rts);
        
        for(int n_scales = scales.getMin(); n_scales <= scales.getMax(); ++n_scales) {
            
            auto part = runForScale(n_response_channels,
                                    n_audio_channels,
                                    n_audio_frames_per_cb,
                                    zero_latency_response_size,
                                    n_scales,
                                    frame_rate,
                                    max_avg_time_per_sample,
                                    os);
            if(!part) {
                os << "Discard n_scales " << n_scales << std::endl;
                continue;
            }
            
            if(!res || (res->getCost() > part->getCost())) {
                res = part;
            }
            
            if(!res) {
                continue;
            }
            if(res->getCost() < max_avg_time_per_sample) {
                os << "Optimization criteria met with " << n_scales << " scaling levels." << std::endl;
                return *res;
            }
            os << "cost " << res->getCost() << " >= " << max_avg_time_per_sample << std::endl;
        }
        throw std::runtime_error("Optimization criteria not met, there are not enough scaling levels.");
        
    }
private:
    static PS runForScale(int n_channels,
                          int const n_audio_channels,
                          int const n_audio_frames_per_cb,
                          int const zero_latency_response_size,
                          int const n_scales,
                          double const frame_rate,
                          double const max_avg_time_per_sample,
                          std::ostream & os)
    {
        auto getEarlyHandlerParams = [&](int countEarlyHandlerCoeffs,
                                         XFFtCostFactors const & xFftFactors) {
            return PartitionAlgo<EarlyHandler>::run(n_channels,
                                                    n_audio_channels,
                                                    n_audio_frames_per_cb,
                                                    countEarlyHandlerCoeffs,
                                                    frame_rate,
                                                    max_avg_time_per_sample,
                                                    os,
                                                    xFftFactors);
            
        };
        if(zero_latency_response_size <= minLatencyLateHandlerWhenEarlyHandlerIsDefaultOptimizedFIRFilter.toInteger()) {
            // in that case there is no late handling at all.
            if(n_scales > 1) {
                LG(WARN, "not enough coefficients to use %d scales", n_scales);
                return {};
            }
            PS o;
            XFFtCostFactors unbiasedXFftCostFactors;
            auto earlyRes = getEarlyHandlerParams(zero_latency_response_size,
                                                  unbiasedXFftCostFactors);
            if(earlyRes) {
                o = mkSubsamplingSetupParams<T, Allocator, FFTTag>(*earlyRes,
                                                                   FinegrainedSetupParam::makeInactive(),
                                                                   1,
                                                                   0);
                o->setCost(earlyRes->getCost());
            }
            return {{o}};
        }
        auto lateRes = PartitionAlgo<EarliestDeepesLateHandler>::run(n_channels,
                                                                     n_scales,
                                                                     n_audio_frames_per_cb,
                                                                     zero_latency_response_size,
                                                                     os);
        if(lateRes) {
            // this fft comes "for free" because the latehandler uses it:
            XFFtCostFactors xFftCostFactors{{
                {lateRes->partition_size, 0.f}
            }};
            // scaleSz might be a little bigger than the scale used during optimization, to have a round number of partitions:
            int const scaleSz = lateRes->partition_size * lateRes->partition_count;
            auto earlyRes = getEarlyHandlerParams(EarliestDeepesLateHandler::nCoefficientsFadeIn +
                                                  EarliestDeepesLateHandler::getLatencyForPartitionSize(lateRes->partition_size).toInteger(),
                                                  xFftCostFactors);
            if(earlyRes) {
                SetupParam ps = mkSubsamplingSetupParams<T, Allocator, FFTTag>(*earlyRes,
                                                                               *lateRes,
                                                                               n_scales,
                                                                               scaleSz);
                ps.setCost(earlyRes->getCost() +
                           lateRes->getCost());
                // we need to adjust because if we have several scales, the last one my contain too many partitions.
                if(ps.handlesCoefficients()) {
                    ps.adjustWork(zero_latency_response_size);
                }
                return ps;
            }
        }
        return {};
    }
};

enum class AudioProcessing {
    Callback, // here, we want 0 latency and the smallest worst case cost per audio callback.
    Offline // here, we want the smallest averaged cost per sample.
};

namespace detail {
template<typename T, AudioProcessing P>
struct OptimalFilter_;

template<typename T>
struct OptimalFilter_<T, AudioProcessing::Callback> {
    // TODO reevaluate once we know when async is better than sync
    using type = ZeroLatencyScaledFineGrainedPartitionnedConvolution<T, a64::Alloc, fft::Fastest>;
};

template<typename T>
struct OptimalFilter_<T, AudioProcessing::Offline> {
    using type = OptimizedFIRFilter<T, a64::Alloc, fft::Fastest>;
};

}

template<typename T, AudioProcessing P>
using ZeroLatencyFilter = typename detail::OptimalFilter_<T,P>::type;

static inline std::ostream & operator << (std::ostream & o, SimulationPhasing const & p)
{
    if(p.mode == SimulationPhasingMode::Off) {
        o << "off";
    }
    else {
        o << "with group size : " << *p.groupSize;
    }
    return o;
}

}
