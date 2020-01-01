
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

  template<typename T, typename FFTTag>
  using OptimizedFIRFilter =
    SplitConvolution <
  // handle the first coefficients using brute force convolution (memory locality makes it faster than anything else):
      FIRFilter<T>,
  // handle subsequent coefficients using FFTs:
      CustomScaleConvolution<FFTConvolutionIntermediate < PartitionnedFFTConvolutionCRTP<T, FFTTag> >>
    >;
  
  template<typename T, typename FFTTag, LatencySemantic Lat = LatencySemantic::DiracPeak>
  using ZeroLatencyScaledFineGrainedPartitionnedSubsampledConvolution =
    SplitConvolution<
  // handle the first coefficients using a zero-latency filter:
      OptimizedFIRFilter<T, FFTTag>,
  // handle subsequent coefficients using a filter optimized for worst audio callback cost:
  SplitConvolution<
    FinegrainedPartitionnedFFTConvolution<T, FFTTag>, // full resolution
    SubSampled<
      Lat,
      Delayed<

  SplitConvolution<
    FinegrainedPartitionnedFFTConvolution<T, FFTTag>, // half resolution
    SubSampled<
      Lat,
      Delayed<

  SplitConvolution<
    FinegrainedPartitionnedFFTConvolution<T, FFTTag>, // quarter resolution
    SubSampled<
      Lat,
      Delayed<

  FinegrainedPartitionnedFFTConvolution<T, FFTTag> // heighth resolution
  
      >
    >
  >
  
      >
    >
  >

      >
    >
  >

  >;

      template<typename T, typename FFTTag>
      using ZeroLatencyScaledFineGrainedPartitionnedConvolution =
        SplitConvolution<
      // handle the first coefficients using a zero-latency filter:
          OptimizedFIRFilter<T, FFTTag>,
      // handle subsequent coefficients using a filter optimized for worst audio callback cost:
          FinegrainedPartitionnedFFTConvolution<T, FFTTag>
      >;
      template<typename T, typename FFTTag, PolicyOnWorkerTooSlow OnWorkerTooSlow>
      using ZeroLatencyScaledAsyncConvolution =
      SplitConvolution<
        OptimizedFIRFilter<T, FFTTag>,
        AsyncCPUConvolution <
            CustomScaleConvolution<FFTConvolutionIntermediate < PartitionnedFFTConvolutionCRTP<T, FFTTag> >>,
            OnWorkerTooSlow
        >
      >;
      template<typename T, typename FFTTag, PolicyOnWorkerTooSlow OnWorkerTooSlow>
      using ZeroLatencyScaledAsyncConvolutionOptimized =
      SplitConvolution<
        ZeroLatencyScaledFineGrainedPartitionnedConvolution<T, FFTTag>,
        AsyncCPUConvolution <
            CustomScaleConvolution<FFTConvolutionIntermediate < PartitionnedFFTConvolutionCRTP<T, FFTTag> >>,
            OnWorkerTooSlow
        >
      >;

template<typename A>
struct PartitionAlgo<CustomScaleConvolution<A>> {
    using Convolution = CustomScaleConvolution<A>;
    using SetupParam = typename Convolution::SetupParam;
    using ScalingParam = typename SetupParam::ScalingParam;
    using PS = PartitionningSpec<SetupParam>;
    using PSpecs = PartitionningSpecs<SetupParam>;

    static PSpecs run(int const n_channels,
                      int const n_audio_channels,
                      int const n_audio_frames_per_cb,
                      int const total_response_size,
                      double const frame_rate,
                      double const max_avg_time_per_sample,
                      std::ostream & os,
                      int const firstSz)
    {
        std::optional<std::pair<std::vector<Scaling>, double>> best =
        getOptimalScalingScheme_ForTotalCpu_ByVirtualSimulation<Convolution>(firstSz, total_response_size);
        
        PS minimalPs;
        if(best) {
            auto scalingParams = scalingsToParams<ScalingParam>(best->first);
            
            minimalPs.cost = SetupParam{scalingParams};
            minimalPs.cost->setCost(best->second);
        }
        return {
            minimalPs,
            minimalPs
        };
    }
};

      template<typename T, typename FFTTag>
      struct PartitionAlgo< OptimizedFIRFilter<T, FFTTag> > {
          using Convolution = OptimizedFIRFilter<T, FFTTag>;
          using LateHandler = typename Convolution::LateHandler;
          using EarlyHandler = typename Convolution::EarlyHandler;
          using SetupParam = typename Convolution::SetupParam;
          using PS = PartitionningSpec<SetupParam>;
          using PSpecs = PartitionningSpecs<SetupParam>;
          
          static PSpecs run(int n_channels,
                            int n_audio_channels,
                            int n_audio_frames_per_cb,
                            int total_response_size,
                            double frame_rate,
                            double max_avg_time_per_sample,
                            std::ostream & os)
          {
              auto const nAsyncScalesDropped = ScaleConvolution_::nDroppedOptimalFor_Split_Bruteforce_Fft;

              int const nEarlyCoeffs = std::min(total_response_size,
                                                static_cast<int>(pow2(nAsyncScalesDropped.toInteger())-1));
              int const nLateCoeffs = std::max(0,
                                               total_response_size - nEarlyCoeffs);
              int const firstSz = static_cast<int>(pow2(nAsyncScalesDropped.toInteger()));

              auto indent = std::make_unique<IndentingOStreambuf>(os);
              auto pSpecsLate = PartitionAlgo<LateHandler>::run(n_channels,
                                                                n_audio_channels,
                                                                n_audio_frames_per_cb,
                                                                nLateCoeffs,
                                                                frame_rate,
                                                                max_avg_time_per_sample,
                                                                os,
                                                                firstSz
                                                                ).getWithSpread();
              indent.reset();
              
              PS ps;
              if(pSpecsLate.cost) {
                  ps.cost = {{}, *pSpecsLate.cost};
                  ps.cost->setCost(pSpecsLate.getCost());
              }
              return {ps,{}};
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
      
      static inline std::ostream & operator << (std::ostream & o, SimulationPhasing const & p) {
          if(p.mode == SimulationPhasingMode::Off) {
            o << "off";
          }
          else {
            o << "with group size : " << *p.groupSize;
          }
        return o;
      }

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
            using SetupParam = typename Convolution::SetupParam;
            using PS = PartitionningSpec<SetupParam>;
            using PSpecs = PartitionningSpecs<SetupParam>;

            // to account for a system that is under heavier load
            static constexpr int safeFactor = 4;

            static PSpecs run(int n_channels,
                              int n_audio_channels,
                              int n_audio_frames_per_cb,
                              int total_response_size,
                              double frame_rate,
                              double max_avg_time_per_sample,
                              std::ostream & os,
                              SimulationPhasing const & phasing,
                              CountDroppedScales const & nAsyncScalesDropped)
            {
                if(n_channels <= 0) {
                    throw std::runtime_error("0 channel");
                }
                if(n_audio_frames_per_cb <= 0) {
                    throw std::runtime_error("0 n_audio_frames_per_cb");
                }
                int const submission_period = n_audio_frames_per_cb;
                int const nMinDroppedCoeffs = n_audio_frames_per_cb;
                int const nLateCoeffs = total_response_size-nMinDroppedCoeffs;
                std::optional<int> maybeQueueSz;
                std::optional<std::pair<std::vector<Scaling>, double>> optimalScaling;
                if(nLateCoeffs <= 0) {
                    // the size of the queue doesn't really matter since there is no late coefficient.
                    maybeQueueSz = 0;
                }
                else {
                    // todo il faudrait peut-etre revoir cela a la hausse
                    // (quitte a augmenter la latence du late) pour reduire la somme total de calculs.
                    int firstSz = static_cast<int>(pow2(nAsyncScalesDropped.toInteger()));
              
                    // Because of the queueing mechanism,
                    // the real number of late coeffs will be _less_ than what this is optimal for.
                    // Hence, as soon as we know the size of the queue, we will be able to deduce the exact
                    // number of late coefficients and adapt the scaling scheme accordingly.
                    //
                    // This future adaptation must be minor : removing work is ok but changing the structure is not ok
                    // because we would have to run the simulation again, as the peak computation times could be worse.
                    optimalScaling = getOptimalScalingScheme_ForTotalCpu_ByVirtualSimulation<AsyncPart>(firstSz, nLateCoeffs);
                    if(!optimalScaling) {
                        throw std::runtime_error("no best");
                    }

                    maybeQueueSz = computeQueueMetricsWithVirtualSimulation(n_channels,
                                                                            n_audio_channels,
                                                                            n_audio_frames_per_cb,
                                                                            nLateCoeffs,
                                                                            frame_rate,
                                                                            optimalScaling->first,
                                                                            submission_period,
                                                                            phasing,
                                                                            os);
                }

                if(!maybeQueueSz) {
                    os << "Synchronous thread is too slow" << std::endl;
                    return {};
                }
                // queueSize == 0 means 0 late coefficient
                int const queueSize = *maybeQueueSz;
                os << "queueSize=" << queueSize << std::endl;
                
                using ScalingParam = typename AsyncPart::SetupParam::ScalingParam;
                  auto scalingParams = scalingsToParams<ScalingParam>([&optimalScaling]() -> std::vector<Scaling> {
                      if(optimalScaling) {
                          return optimalScaling->first;
                      }
                      else {
                          return {};
                      }
                  }());

                PS minimalPs;
                minimalPs.cost = SetupParam{
                    submission_period,
                    queueSize,
                    {
                        scalingParams
                    }
                };
                minimalPs.cost->setCost(optimalScaling ? optimalScaling->second : 0.);
                return {
                    minimalPs,
                    minimalPs
                };
            }
        private:
            
            static std::optional<int>
            computeQueueMetricsWithVirtualSimulation(int const n_channels,
                                                     int const n_audio_channels,
                                                     int const n_audio_frames_per_cb,
                                                     int const nLateCoeffs,
                                                     double const frame_rate,
                                                     std::vector<Scaling> const & scaling,
                                                     int const submission_period_frames,
                                                     SimulationPhasing const & phasing,
                                                     std::ostream & os)
            {
                using Sim = typename AsyncPart::Simulation;
                
                std::vector<std::unique_ptr<Sim>> simu_convs;
                simu_convs.reserve(n_channels);
                for(int i=0; i<n_channels; ++i) {
                    simu_convs.push_back(std::make_unique<Sim>());
                }

                PS minimalPs;
                using ScalingParam = typename Sim::SetupParam::ScalingParam;
                auto scalingParams = scalingsToParams<ScalingParam>(scaling);
              
                for(auto & c : simu_convs) {
                    c->setup({scalingParams});
                    c->setCoefficientsCount(nLateCoeffs);
                }
                if(phasing.mode == SimulationPhasingMode::On)
                {
                    int n=0;
                    int const total_sz = phasing.groupSize ? *phasing.groupSize : simu_convs.size();
                    Assert(total_sz);
                    typename Sim::SetupParam sp(
                        scalingParams // not used
                    );
                    for(auto & c : simu_convs) {
                        dephase(total_sz,
                                n % total_sz,
                                sp,
                                *c);
                        ++n;
                    }
                }
                
                // at the moment, every AsyncCPUConvolution has its own worker thread, so we simulate
                // 'n_channels' threads runing at the same time, to see what percentage of the time they
                // are actually running on the cpu, and what their max 'sleep' time is.
                
                auto const nThreads = n_channels;
                os << "With " << nThreads << " threads:" << std::endl;
                auto indent = std::make_unique<IndentingOStreambuf>(os);

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

                int queueSize = computeQueueSize([&simu_convs,
                                                  min_allocation_factor](){
                    // take only the first one into account (each of them is equivalent over a period)
                    return simu_convs[0]->simuStep() / min_allocation_factor;
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
                    SetupParam sp{
                        1,
                        1,
                        {
                            CountDroppedScales(1)
                        }
                    };
                    for(auto & c : async_convs) {
                        dephase(total_sz,
                                n % total_sz,
                                sp,
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
      
      template<typename T, typename FFTTag, PolicyOnWorkerTooSlow OnWorkerTooSlow>
      struct PartitionAlgo< ZeroLatencyScaledAsyncConvolution<T, FFTTag, OnWorkerTooSlow> > {
          using Convolution = ZeroLatencyScaledAsyncConvolution<T, FFTTag, OnWorkerTooSlow>;
          using SetupParam = typename Convolution::SetupParam;
          using PS = PartitionningSpec<SetupParam>;
          using PSpecs = PartitionningSpecs<SetupParam>;

          using LateHandler = typename Convolution::LateHandler;
          using EarlyHandler = typename Convolution::EarlyHandler;

          static PSpecs run(int n_channels,
                            int n_audio_channels,
                            int n_audio_frames_per_cb,
                            int total_response_size,
                            double frame_rate,
                            double max_avg_time_per_sample,
                            std::ostream & os,
                            SimulationPhasing const & phasing)
        {
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
                                                              total_response_size,
                                                              frame_rate,
                                                              max_avg_time_per_sample,
                                                              os,
                                                              phasing,
                                                              nAsyncScalesDropped).getWithSpread();
            PS ps;
            if(pSpecsLate.cost) {
                auto nEarlyHandlerCoeffs = pSpecsLate.cost->getImpliedLatency();
                auto pSpecsEarly = PartitionAlgo<EarlyHandler>::run(n_channels,
                                                                  n_audio_channels,
                                                                  n_audio_frames_per_cb,
                                                                  nEarlyHandlerCoeffs,
                                                                  frame_rate,
                                                                  max_avg_time_per_sample,
                                                                  os).getWithSpread();
                if(pSpecsEarly.cost) {
                    ps.cost = {
                        *pSpecsEarly.cost,
                        *pSpecsLate.cost
                    };
                    ps.cost->setCost(pSpecsEarly.getCost() +
                                     pSpecsLate.getCost());
                }
            }
            return {ps,{}};
        }
      };
      template<typename T, typename FFTTag, PolicyOnWorkerTooSlow OnWorkerTooSlow>
      struct PartitionAlgo< ZeroLatencyScaledAsyncConvolutionOptimized<T, FFTTag, OnWorkerTooSlow> > {
          using Convolution = ZeroLatencyScaledAsyncConvolutionOptimized<T, FFTTag, OnWorkerTooSlow>;
          using SetupParam = typename Convolution::SetupParam;
          using PS = PartitionningSpec<SetupParam>;
          using PSpecs = PartitionningSpecs<SetupParam>;

          using EarlyHandler = typename Convolution::EarlyHandler;
          using LateHandler = typename Convolution::LateHandler;

          using LateSetupParam = typename LateHandler::SetupParam;
      
      static PSpecs run(int n_channels,
                        int n_audio_channels,
                        int n_audio_frames_per_cb,
                        int total_response_size,
                        double frame_rate,
                        double max_avg_time_per_sample,
                        std::ostream & os,
                        SimulationPhasing const & phasing)
      {
          os << "Late handler optimization:" << std::endl;

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

          auto indent = std::make_unique<IndentingOStreambuf>(os);
          auto pSpecsLate = PartitionAlgo<LateHandler>::run(n_channels,
                                                            n_audio_channels,
                                                            n_audio_frames_per_cb,
                                                            total_response_size,
                                                            frame_rate,
                                                            max_avg_time_per_sample,
                                                            os,
                                                            phasing,
                                                            nAsyncScalesDropped).getWithSpread();
          indent.reset();
          
          PS ps;
          if(pSpecsLate.cost) {
              int nEarlyCoeffs = deduceEarlyHandlerCoeffs(*(pSpecsLate.cost));
              os << "Deduced early handler coeffs:" << nEarlyCoeffs << std::endl;
              os << "Early handler optimization:" << std::endl;
              
              auto indent = std::make_unique<IndentingOStreambuf>(os);
              auto pSpecsEarly = PartitionAlgo<EarlyHandler>::run(n_channels,
                                                                  n_audio_channels,
                                                                  n_audio_frames_per_cb,
                                                                  nEarlyCoeffs,
                                                                  frame_rate,
                                                                  max_avg_time_per_sample,
                                                                  os).getWithSpread();
              indent.reset();
              
              if(pSpecsEarly.cost) {
                  ps.cost = {*pSpecsEarly.cost, *pSpecsLate.cost};
                  ps.cost->setCost(pSpecsEarly.getCost() + pSpecsLate.getCost());
              }
          }
        return {ps,{}};
      }
      private:
      static int deduceEarlyHandlerCoeffs(LateSetupParam p) {
        return p.getImpliedLatency();
      }
      };

namespace SameSizeScales {
    /*
     Assuming that scales have resolutions of 1,2,4,8...
     for every scale, the number of blocks = { N(k), k in 1 .. S }
     Assuming that between scale i and (i+1) there are nOverlap*2^(i-1) blocks in common
     N = sum (n in 0 .. (S-1)) 2^n N(n+1) - scaleFadeSz::inSmallerUnits * sum (n in 1.. (S-1)) 2^(n-1)
     and since we want the number of blocks to be equal
     (hence the namespace name 'SameSizeScales'), we will solve this:
     N = NF * sum (n in 0 .. (S-1)) 2^n - scaleFadeSz::inSmallerUnits * sum (n in 0.. (S-2)) 2^n
     N = NF * (2^S - 1) - scaleFadeSz::inSmallerUnits * (2^(S-1)-1)
     NF = (N + scaleFadeSz::inSmallerUnits * (2^(S-1)-1)) / (2^S - 1)
     */
    static inline int get_scale_sz(int response_sz, int n_scales) {
      assert(response_sz>=0);
      int numerator = response_sz + static_cast<int>(scaleFadeSz::inSmallerUnits * (static_cast<int>(pow2(n_scales-1)) - 1));
      int denominator = static_cast<int>(pow2(n_scales)) - 1;
      int res = ceil(numerator / static_cast<double>(denominator));
      if(subSamplingAllowsEvenNumberOfCoefficients || n_scales <= 1) {
        return res;
      }
      if(res%2) {
        return res + 1;
      }
      return res;
    }
    static inline int get_max_response_sz(int n_scales, int sz_one_scale) {
      return
      sz_one_scale * (pow2(n_scales) - 1) -
      scaleFadeSz::inSmallerUnits * (static_cast<int>(pow2(n_scales-1)) - 1);
    }
    
    int constexpr getDelays(int scale_sz, int partition_sz) {
      // split = nFadeIn + latB - latA
      // scale_sz = nFadeIn + (1 + 2*(delay + latA)) - latA
      // scale_sz - 1 - nFadeIn = 2*delay + latA
      // delay = 0.5 * (scale_sz - 1 - nFadeIn - latA)
      // delay = 0.5 * (scale_sz - 1 - nFadeIn - (2*partition_sz-1))
      // delay = (0.5 * (scale_sz - nFadeIn)) - partition_sz
      assert(0 == (scale_sz-scaleFadeSz::inSmallerUnits) % 2);
      return ((scale_sz-scaleFadeSz::inSmallerUnits) / 2) - partition_sz;
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
  template<typename A>
  struct EarlyestDeepest< Delayed<A> > {
    using type = typename EarlyestDeepest<A>::type;
  };
  template<LatencySemantic Lat, typename A>
  struct EarlyestDeepest< SubSampled<Lat, A> > {
    using type = typename EarlyestDeepest<A>::type;
  };

  
  constexpr int minLatencyLateHandlerWhenEarlyHandlerIsDefaultOptimizedFIRFilter = ScaleConvolution_::latencyForDroppedConvolutions(ScaleConvolution_::nDroppedOptimalFor_Split_Bruteforce_Fft);

  template<typename C>
  int constexpr lateHandlerLatency(int partition_sz) {
    using LateHandler = typename EarlyestDeepest<typename C::LateHandler>::type;
    auto res = LateHandler::getLatencyForPartitionSize(partition_sz);
    return res;
  }
  
  template<typename C>
  int constexpr getLateHandlerMinLg2PartitionSz() {
    using LateHandler = typename EarlyestDeepest<typename C::LateHandler>::type;

    int partition_sz = 1;
    for(;;partition_sz *= 2) {
      if(LateHandler::getLatencyForPartitionSize(partition_sz) >= minLatencyLateHandlerWhenEarlyHandlerIsDefaultOptimizedFIRFilter) {
        break;
      }
    }
    
    return power_of_two_exponent(partition_sz);
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
      

      template<typename T, typename FFTTag, typename SPEarly, typename SPLate>
      auto mkSubsamplingSetupParams(SPEarly const & early_params,
                                    SPLate const & late_params,
                                    int const n_scales,
                                    int const scale_sz) {
        using SetupParam = typename ZeroLatencyScaledFineGrainedPartitionnedSubsampledConvolution<T,FFTTag>::SetupParam;
        int delay = 0;
        if(n_scales > 1) {
          delay = SameSizeScales::getDelays(scale_sz, late_params.partition_size);
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
              (n_scales >= 2)?delay:0,
              {
                ps[1],
                {
                  (n_scales >= 3)?delay:0,
                  {
                    ps[2],
                    {
                      (n_scales >= 4)?delay:0,
                      ps[3]
                    }
                  }
                }
              }
            }
          }
        };
      }

      
  template<typename T, typename FFTTag>
  struct PartitionAlgo< ZeroLatencyScaledFineGrainedPartitionnedConvolution<T,FFTTag> > {
    using NonAtomicConvolution = ZeroLatencyScaledFineGrainedPartitionnedConvolution<T,FFTTag>;
    
  private:
    using LateHandler = typename NonAtomicConvolution::LateHandler;
    using EarlyHandler = typename NonAtomicConvolution::EarlyHandler;
  public:
      using SetupParam = typename NonAtomicConvolution::SetupParam;
    using PS = PartitionningSpec<SetupParam>;
    using PSpecs = PartitionningSpecs<SetupParam>;

    static PSpecs run(int n_channels,
                      int n_audio_channels,
                      int n_audio_frames_per_cb,
                      int total_response_size,
                      double frame_rate,
                      double max_avg_time_per_sample,
                      std::ostream & os) {
      auto getEarlyHandlerParams = [&](int countEarlyHandlerCoeffs) {
          return PartitionAlgo<EarlyHandler>::run(n_channels,
                                                  n_audio_channels,
                                                  n_audio_frames_per_cb,
                                                  countEarlyHandlerCoeffs,
                                                  frame_rate,
                                                  max_avg_time_per_sample,
                                                  os).getWithSpread();
          
      };
      if(total_response_size <= minLatencyLateHandlerWhenEarlyHandlerIsDefaultOptimizedFIRFilter) {
        // in that case there is no late handling at all.
        PS ps;
        auto earlyRes = getEarlyHandlerParams(total_response_size);
        if(earlyRes.cost)
        {
          ps.cost =
          {
            *earlyRes.cost,
            FinegrainedSetupParam::makeInactive()
          };
          ps.cost->setCost(earlyRes.cost->getCost());
        }
        return {ps,ps};
      }
      auto late_handler_response_size_for_partition_size =
      [total_response_size](int partition_size) -> Optional<int> {
        if(LateHandler::getLatencyForPartitionSize(partition_size) < minLatencyLateHandlerWhenEarlyHandlerIsDefaultOptimizedFIRFilter) {
          // invalid case. We pass 'getMinLg2PartitionSz()' to 'run' so that
          // this case doesn't happen on the first try.
          return {};
        }
        // we substract the count of coefficients handled by the early coefficients handler.
        // it favors long partitions because we don't take into account
        // the cost of the early coefficient handler, but for long responses, where optimization matters,
        // the induced bias is negligible.
        int n_coeffs_early_handler = lateHandlerLatency<NonAtomicConvolution>(partition_size);
        assert(n_coeffs_early_handler >= minLatencyLateHandlerWhenEarlyHandlerIsDefaultOptimizedFIRFilter);

        auto late_response_sz = std::max(0,total_response_size - n_coeffs_early_handler);

        if(late_response_sz <= 0) {
      // should we return {0} ?
          return {};
        }
        return {late_response_sz};
      };
      auto lateRes = PartitionAlgo<LateHandler>::run(n_channels,
                                                     1,
                                                     n_audio_frames_per_cb,
                                                     late_handler_response_size_for_partition_size,
                                                     getLateHandlerMinLg2PartitionSz<NonAtomicConvolution>(),
                                                     os);
      PS withSpread;
      if(lateRes.with_spread.cost) {
          auto earlyRes = getEarlyHandlerParams(LateHandler::nCoefficientsFadeIn +
                                                LateHandler::getLatencyForPartitionSize(lateRes.with_spread.cost.value().partition_size));
          if(earlyRes.cost) {
              withSpread.cost = SetupParam
              {
                  *earlyRes.cost,
                  lateRes.with_spread.cost.value()
              };
              withSpread.cost->setCost(earlyRes.cost->getCost() +
                                       lateRes.with_spread.cost->getCost());
          }
      }
      PS withoutSpread;
      if(lateRes.without_spread.cost) {
          auto earlyRes = getEarlyHandlerParams(LateHandler::nCoefficientsFadeIn +
                                                LateHandler::getLatencyForPartitionSize(lateRes.without_spread.cost.value().partition_size));
          if(earlyRes.cost) {
              withoutSpread.cost = SetupParam
              {
                  *earlyRes.cost,
                  lateRes.without_spread.cost.value()
              };
              withoutSpread.cost->setCost(earlyRes.cost->getCost() +
                                          lateRes.without_spread.cost->getCost());
          }
      }
      return {
        withSpread,
        withoutSpread
      };
    }
  };
  
      template<typename T>
      struct WithNScales {
        T val;
        int n_scales;
      
        auto getCost() const { return val.getCost(); }
      
      void logReport(int n_channels,
                     double theoretical_max_avg_time_per_frame,
                     std::ostream & os) {
          val.logReport(n_channels, theoretical_max_avg_time_per_frame, os);
          os << "- using ";
          // TODO range of subsampling regions: highest quality region: xxx samples / 2-subsampled region : xxx samples / 4-subsampled
          if(n_scales <= 1) {
            os << "full tail resolution";
          }
          else {
            os << "reduced tail resolution with " << n_scales - 1 << " subsampling regions";
          }
      }
      };
      
      template<typename T, typename FFTTag>
      struct PartitionAlgo< ZeroLatencyScaledFineGrainedPartitionnedSubsampledConvolution<T,FFTTag> > {
        using NonAtomicConvolution = ZeroLatencyScaledFineGrainedPartitionnedSubsampledConvolution<T,FFTTag>;
        
      private:
        using EarlyHandler = typename NonAtomicConvolution::EarlyHandler;
        using LateHandler = typename EarlyestDeepest<typename NonAtomicConvolution::LateHandler>::type;
      public:
        using SetupParam = typename NonAtomicConvolution::SetupParam;//LateHandler::SetupParam;
        using PS = PartitionningSpec<SetupParam>;
        using PSpecs = PartitionningSpecs2<WithNScales<SetupParam>>;

      // il faudrait que a chaque fois ca renvoie PartitionningSpecs<NonAtomicConvolution::SetupParam> pour bien generaliser
        static PSpecs run(int const n_response_channels,
                          int const n_audio_channels,
                          int const n_audio_frames_per_cb,
                          int const n_response_frames,
                          double const frame_rate,
                          double const max_avg_time_per_sample,
                          std::ostream & os,
                          ResponseTailSubsampling rts) {
      
      Optional<PSpecs> res;
      
      range<int> const scales = getScaleCountRanges<NonAtomicConvolution>(rts);

      for(int n_scales = scales.getMin(); n_scales <= scales.getMax(); ++n_scales) {
        
        auto partit = runForScale(n_response_channels,
                                 n_audio_channels,
                                 n_audio_frames_per_cb,
                                 n_response_frames,
                                 n_scales,
                                 frame_rate,
                                  max_avg_time_per_sample,
                                 os);
        auto & part = partit.getWithSpread();
        if(!part.cost) {
          os << "Discard n_scales " << n_scales << std::endl;
          continue;
        }
        
        if(!res || (res->getCost() && (*res->getCost() > part.getCost()))) {
          res = {part, part};
        }

        if(!res || !res->getCost()) {
            continue;
        }
        if(*res->getCost() < max_avg_time_per_sample) {
          os << "Optimization criteria met with " << n_scales << " scaling levels." << std::endl;
          return *res;
        }
        os << "cost " << *res->getCost() << " >= " << max_avg_time_per_sample << std::endl;
      }
      throw std::runtime_error("Optimization criteria not met, there are not enough scaling levels.");

      }
  private:
      static PSpecs runForScale(int n_channels,
                   int const n_audio_channels,
                   int const n_audio_frames_per_cb,
                   int const total_response_size,
                   int const n_scales,
                   double const frame_rate,
                   double const max_avg_time_per_sample,
                   std::ostream & os)
      {
          auto getEarlyHandlerParams = [&](int countEarlyHandlerCoeffs) {
              return PartitionAlgo<EarlyHandler>::run(n_channels,
                                                      n_audio_channels,
                                                      n_audio_frames_per_cb,
                                                      countEarlyHandlerCoeffs,
                                                      frame_rate,
                                                      max_avg_time_per_sample,
                                                      os).getWithSpread();
              
          };
          if(total_response_size <= minLatencyLateHandlerWhenEarlyHandlerIsDefaultOptimizedFIRFilter) {
            // in that case there is no late handling at all.
            if(n_scales > 1) {
              LG(WARN, "not enough coefficients to use %d scales", n_scales);
              return {};
            }
            std::optional<WithNScales<SetupParam>> o;
            auto earlyRes = getEarlyHandlerParams(total_response_size);
            if(earlyRes.cost) {
              o = {
                mkSubsamplingSetupParams<T, FFTTag>(*earlyRes.cost,
                                                    FinegrainedSetupParam::makeInactive(),
                                                    1,
                                                    0),
                n_scales
              };
              o->val.setCost(earlyRes.cost->getCost());
            }
            return {{o},{o}};
          }
          auto scale_size_for_partition_size =
          [total_response_size, n_scales](int partition_size) -> Optional<int> {
            if(LateHandler::getLatencyForPartitionSize(partition_size) < minLatencyLateHandlerWhenEarlyHandlerIsDefaultOptimizedFIRFilter) {
              // invalid case. We pass 'getMinLg2PartitionSz()' to 'run' so that
              // this case doesn't happen on the first try.
              return {};
            }
            // we substract the count of coefficients handled by the early coefficients handler.
            // it favors long partitions because we don't take into account
            // the cost of the early coefficient handler, but for long responses, where optimization matters,
            // the induced bias is negligible.
            int n_coeffs_early_handler = lateHandlerLatency<NonAtomicConvolution>(partition_size);
            assert(n_coeffs_early_handler >= minLatencyLateHandlerWhenEarlyHandlerIsDefaultOptimizedFIRFilter);

            auto late_response_sz = std::max(0,total_response_size - n_coeffs_early_handler);

            auto scale_sz = SameSizeScales::get_scale_sz(late_response_sz, n_scales);
            if(scale_sz <= 0) {
              return {};
            }
            // TODO return also the size of early coefficients, and take the early cost into account.
            if(n_scales > 1) {
              // verify that we can have more than one scale (i.e the delay needed to scale is strictly positive)
              // in theory we could relax the constraint (0 delays are ok but the implementation
              // doesn't support that).
              if(SameSizeScales::getDelays(scale_sz, partition_size) <= 0) {
                return {};
              }
            }
            return {scale_sz};
          };
          auto lateRes = PartitionAlgo<LateHandler>::run(n_channels,
                                                         n_scales,
                                                         n_audio_frames_per_cb,
                                                         scale_size_for_partition_size,
                                                         getLateHandlerMinLg2PartitionSz<NonAtomicConvolution>(),
                                                         os);
          std::optional<WithNScales<SetupParam>> withSpread;
          if(lateRes.with_spread.cost) {
            std::optional<int> const mayScaleSz = scale_size_for_partition_size(lateRes.with_spread.cost->partition_size);
            if(mayScaleSz) {
              auto earlyRes = getEarlyHandlerParams(LateHandler::nCoefficientsFadeIn +
                                                    LateHandler::getLatencyForPartitionSize(lateRes.with_spread.cost.value().partition_size));
              if(earlyRes.cost) {
                  withSpread = {
                      mkSubsamplingSetupParams<T, FFTTag>(*earlyRes.cost,
                                                          lateRes.with_spread.cost.value(),
                                                          n_scales,
                                                          mayScaleSz.value()),
                      n_scales
                  };
                  withSpread->val.setCost(lateRes.with_spread.cost->getCost() +
                                          earlyRes.cost->getCost());
              }
            }
            else {
                throw std::logic_error("no scale size");
            }
          }
          std::optional<WithNScales<SetupParam>> withoutSpread;
          if(lateRes.without_spread.cost) {
            std::optional<int> const mayScaleSz = scale_size_for_partition_size(lateRes.without_spread.cost->partition_size);
            if(mayScaleSz) {
              auto earlyRes = getEarlyHandlerParams(LateHandler::nCoefficientsFadeIn +
                                                    LateHandler::getLatencyForPartitionSize(lateRes.without_spread.cost.value().partition_size));
              if(earlyRes.cost) {
                withoutSpread = {
                  mkSubsamplingSetupParams<T, FFTTag>(*earlyRes.cost,
                                                      lateRes.without_spread.cost.value(),
                                                      n_scales,
                                                      mayScaleSz.value()),
                  n_scales
                };
                withoutSpread->val.setCost(lateRes.without_spread.cost->getCost() +
                                           earlyRes.cost->getCost());
              }
            }
          }
          return {
              {withSpread},
              {withoutSpread}
          };
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
      using type = ZeroLatencyScaledFineGrainedPartitionnedConvolution<T, fft::Fastest>;
    };
    
    template<typename T>
    struct OptimalFilter_<T, AudioProcessing::Offline> {
      using type = OptimizedFIRFilter<T, fft::Fastest>;
    };
  }
  
  template<typename T, AudioProcessing P>
  using ZeroLatencyFilter = typename detail::OptimalFilter_<T,P>::type;
 

  template<typename Algo>
  void debugSteps(Algo & c, int nMax=20) {
    using namespace std;

    cout << "--" << endl;

    auto v = c.debugStep(1);
    int k = 0;
    cout << to_string(k) << " : " << to_string(v.first)<< " " << to_string(v.second) << endl;
    ++k;
    for(; k<nMax; ++k) {
      v = c.debugStep(0);
      cout << to_string(k) << " : " << to_string(v.first)<< " " << to_string(v.second) << endl;
    }
  }
}
