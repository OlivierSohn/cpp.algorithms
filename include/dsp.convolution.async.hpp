
namespace imajuscule::audio
{
enum class PolicyOnWorkerTooSlow {
    PermanentlySwitchToDry,
    Wait // for testing purposes
};

struct AsyncCPUConvolutionConstants {
    // We leave room for one more element in worker_2_rt queue so that the worker
    //   can push immediately if its computation is faster than the real-time thread
    //   (which is very unlikely, though).
    static constexpr int queue_room_sz = 1;
};

template<typename InnerParams>
struct AsyncSetupParam : public Cost {
    static constexpr bool has_subsampling = InnerParams::has_subsampling;

    using AsyncParam = InnerParams;
    static constexpr int queue_room_sz = AsyncCPUConvolutionConstants::queue_room_sz;

    AsyncSetupParam(int inputSubmissionPeriod,
               int queueSize,
               InnerParams const & asyncParams)
    : inputSubmissionPeriod(inputSubmissionPeriod)
    , queueSize(queueSize)
    , asyncParams(asyncParams)
    {}
    
    int inputSubmissionPeriod;
    int queueSize;
    InnerParams asyncParams;
    
    bool handlesCoefficients() const {
        return queueSize > 0 && inputSubmissionPeriod > 0 && asyncParams.handlesCoefficients();
    }

    Latency getImpliedLatency() const {
        Assert(handlesCoefficients());
        return asyncParams.getImpliedLatency() +
        Latency( inputSubmissionPeriod*((queueSize-queue_room_sz) + 1) - 1);
    }

    int getBiggestScale() const {
        return inputSubmissionPeriod;
    }

    void adjustWork(int targetNCoeffs) {
        asyncParams.adjustWork(targetNCoeffs);
        if(!asyncParams.handlesCoefficients()) {
            queueSize = 0; // not sure if that is needed
        }
        if(!handlesCoefficients()) {
            setCost(0.);
        }
    }
    
    template<Overlap Mode, typename FFTAlgo>
    MinSizeRequirement getMinSizeRequirement(int const maxVectorSz) const
    {
        // we write y by chunks of submissionPeriod :
        int const y_size = inputSubmissionPeriod * countPartitions(maxVectorSz,
                                                                   inputSubmissionPeriod);

        return {
            maxVectorSz, // x size (pour copier les x dans la file d'attente).
            y_size, // y size
            {}, // ffts : besoin de rien pour l'instant puisque les fft sont recalculées dans la partie async
            0 // work size
        };
    }
    
    void logSubReport(std::ostream & os) const override {
        os << "Async, period : " << inputSubmissionPeriod <<  " size : " << queueSize << std::endl;
        {
            IndentingOStreambuf i(os);
            asyncParams.logSubReport(os);
        }
    }
    
    template<typename F>
    void forEachUsingSameContext(F f) const {
        f(*this);
        // do not recurse into asyncParams because the async part uses another context.
    }
};

}
