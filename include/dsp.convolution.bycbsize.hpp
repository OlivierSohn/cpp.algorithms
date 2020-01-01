

namespace imajuscule {

// todo use std::lerp when c++ 20 is available
template<typename T1, typename T2, typename T3>
auto lerp(T1 from, T2 to, T3 ratio) {
    return from + ratio * (to-from);
}

template<typename Rev>
struct ConvReverbsByBlockSize {
    
    static constexpr int nOut = Rev::nOut;
    static_assert(nOut != 0);
    
    using FPT = typename Rev::FPT;
    
    template<typename ...Args>
    void setConvolutionReverbIR(int const n_sources,
                                DeinterlacedBuffers<FPT> const & deinterlaced,
                                int const maxCountFramesPerBlock,
                                int const minSize,
                                double sampleRate,
                                std::map<int, ConvReverbOptimizationReport> & results,
                                Args... args)
    {        
        // it's ok not to lock here : we are the only thread using these
        reverbsByBlockSize.clear();
        deprecated.clear();
        current = {};
        
        results.clear();

        // we make reverbs optimized exactly for maxCountFramesPerBlock, and smaller powers of 2
        for(int max2 = maxCountFramesPerBlock;
            max2 > 0;
            max2 = floor_power_of_two(max2)/2)
        {
            auto r = std::make_unique<Rev>();

            try {
                ConvReverbOptimizationReportSuccess result;
                std::stringstream os;
                r->setConvolutionReverbIR(n_sources,
                                          deinterlaced,
                                          max2,
                                          sampleRate,
                                          os,
                                          result.structure,
                                          args...);
                result.optimizationReport = os.str();
                {
                    auto res = results.try_emplace(max2, result);
                    assert(res.second);
                }
                {
                    auto res = reverbsByBlockSize.try_emplace(max2, std::move(r));
                    assert(res.second);
                }
            }
            catch (std::exception const & e) {
                auto res = results.try_emplace(max2, e.what());
                assert(res.second);
            }
            
            if(max2 <= minSize) {
                break;
            }
        }
        deprecated.reserve(reverbsByBlockSize.size());
    }

    void abruptlySetConvolutionReverbWetRatio(double wet) {
        transitionConvolutionReverbWetRatio(wet, 0);
    }
    void transitionConvolutionReverbWetRatio(double wet, int duration) {
      wetRatio.smoothly(clamp_ret(wet,0.,1.), duration);
    }
    
      double getWetRatioWithoutStepping() const {
          return wetRatio.getWithoutStepping();
      }
    
    int declareBlockSize(int blocksize) {
        if(!current) {
            current = mkCurrentFromBlockSize(blocksize);
        }
        else if(current->hostBlockSize != blocksize) {
            auto prevCurrent = current;
            current = mkCurrentFromBlockSize(blocksize);
            if(!current) {
                deprecated.clear();
            }
            else if(current->blockSizeHypothesys() != prevCurrent->blockSizeHypothesys()) {
                // if the new current is in deprecated, remove it:
                auto sameKey = [key = current->getKey()](auto const & other){
                    return other.getKey() == key;
                };
                deprecated.erase(std::remove_if(deprecated.begin(), deprecated.end(), sameKey),
                                 deprecated.end());
                
                // add the former current to deprecated
                deprecated.emplace_back(prevCurrent->reverbIt);
            }
        }
        return current ? current->blockSizeHypothesys() : -1;
    }
    

    template<typename FPT2>
    bool apply(FPT2 ** io_buffers,
               int nInputBuffers,
               FPT2 ** work_buffers,
               int nWorkBuffers,
               int nFramesToCompute)
    {
        bool has_step_errors = false;
        
        if(likely(current)) {
            has_step_errors = !current->getReverb().assignWetVectorized(io_buffers,
                                                                       nInputBuffers,
                                                                       work_buffers,
                                                                       nWorkBuffers,
                                                                       nFramesToCompute,
                                                                       nFramesToCompute);
           // by design, deprecated is non empty only if current is non empty
            for(auto & d : deprecated) {
                int const nApply = std::min(d.numFramesToGo, nFramesToCompute);
                d.numFramesToGo -= nApply;
                has_step_errors = !d.getReverb().addWetInputZeroVectorized(work_buffers,
                                                                           nWorkBuffers,
                                                                           nApply,
                                                                           nApply) || has_step_errors;
            }

            // clean up deprecated
            {
                auto isDone = [](auto const & d) {
                    return d.numFramesToGo <= 0;
                };
                deprecated.erase(std::remove_if(deprecated.begin(), deprecated.end(), isDone),
                                 deprecated.end());
            }
        }
        else {
            for(int j=0; j<nWorkBuffers; ++j) {
                fft::RealSignal_<fft::Fastest, FPT2>::zero_n_raw(work_buffers[j], nFramesToCompute);
            }
        }
        
        // apply dry/wet effect
        if(nInputBuffers == nWorkBuffers) {
            for(int i=0; i<nFramesToCompute; ++i) {
                auto wet = wetRatio.step();
                for(int j=0; j<nWorkBuffers; ++j) {
                    auto * wet_buffer = work_buffers[j];
                    auto * dry_buffer = io_buffers[j];
                    dry_buffer[i] = lerp(dry_buffer[i],
                                         wet_buffer[i],
                                         wet);
                }
            }
        }
        
        return has_step_errors;
    }
    
    void flushToSilence() {
        forEachActiveReverb([](auto & rev) {
            rev.flushToSilence();
        });
    }
    
private:
    static constexpr int dryWetPolyOrder = 1; // 1 only because there is no volume change between dry/wet
    using WetRatio = smoothVar<double, dryWetPolyOrder>;
    WetRatio wetRatio = {1.};

    using Map = std::map<int, std::unique_ptr<Rev>>;
    using ConstIterator = typename Map::const_iterator;
    struct Current {
        
        int blockSizeHypothesys() const {
            return reverbIt->first;
        }

        int getKey() const {
            return reverbIt->first;
        }
        Rev const & getReverb() const {
            return *reverbIt->second.get();
        }
        Rev & getReverb() {
            return *reverbIt->second.get();
        }

        int hostBlockSize;
        ConstIterator reverbIt;
    };
    Map reverbsByBlockSize;
    std::optional<Current> current;
    struct DeprecatedReverb {
        DeprecatedReverb(ConstIterator it)
        : reverbIt(it),
        numFramesToGo(it->second->getUnpaddedSize())
        {}
        
        ConstIterator reverbIt;
        int numFramesToGo;

        int getKey() const {
            return reverbIt->first;
        }
        Rev const & getReverb() const {
            return *reverbIt->second.get();
        }
        Rev & getReverb() {
            return *reverbIt->second.get();
        }
    };
    std::vector<DeprecatedReverb> deprecated;


    template<typename F>
    void forEachActiveReverb(F f) {
        if(current) {
            f(current->getReverb());
        }
        for(auto & d : deprecated) {
            f(d.getReverb());
        }
    }
    std::optional<Current> mkCurrentFromBlockSize(int blocksize) const {
        if(reverbsByBlockSize.empty()) {
            return {};
        }
        // conservatively, we use the closest size smaller or equal to blocksize.
        auto it = reverbsByBlockSize.upper_bound(blocksize);
        // now, it points to a reverb that was optimized for a strictly bigger block size
        if(it == reverbsByBlockSize.begin()) {
            return {};
        }
        --it;
        // now, it points to a reverb that was optimized for a smaller or equal block size
        return {{blocksize, it}};
    }
    
};

} // NS imajuscule

