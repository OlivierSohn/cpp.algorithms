

namespace imajuscule::audio {

// todo use std::lerp when c++ 20 is available
template<typename T1, typename T2, typename T3>
auto lerp(T1 from, T2 to, T3 ratio) {
  return from + ratio * (to-from);
}

static inline constexpr int maxVectorSizeFromBlockSizeHypothesis(int const blockSzHypothesis) {
  /*
   When using 'blockSzHypothesis':
   . Callbacks that have a "power of 2" size will use a single vector call.
   . Callbacks that are _not_ a "power of 2" will use 2 vector calls.

   When using 'blockSzHypothesis * 2':
   . All callbacks use a single vector call
   (at the cost of a _slighty_ higher overhead)
   */
  return blockSzHypothesis * 2;
}


template<typename Rev>
struct ConvReverbsByBlockSize {

  static constexpr int nOut = Rev::nOut;
  static_assert(nOut != 0);

  using FPT = typename Rev::FPT;
  using WorkCplxFreqs = typename Rev::WorkCplxFreqs;

  void clear() {
    // it's ok not to lock here : we are the only thread using these
    reverbsByBlockSize.clear();
    memories.clear();
    deprecated.clear();
    current = {};

    reverbUnpaddedLength = 0;
    work = {};
  }

  template<typename ...Args>
  void setConvolutionReverbIR(int const n_sources,
                              DeinterlacedBuffers<FPT> const & deinterlaced,
                              int const maxCountFramesPerBlock,
                              int const minSize,
                              double sampleRate,
                              std::map<int, ConvReverbOptimizationReport> & results,
                              Args... args)
  {
    clear();

    results.clear();
    reverbUnpaddedLength = deinterlaced.countFrames();

    // we make reverbs optimized exactly for maxCountFramesPerBlock, and smaller powers of 2
    for(int max2 = maxCountFramesPerBlock;
        max2 > 0;
        max2 = floor_power_of_two(max2)/2)
    {
      auto memory = std::make_unique<typename Rev::MemResource::type>();

      auto r = std::make_unique<Rev>();

      try {
        ConvReverbOptimizationReportSuccess result;
        std::stringstream os;
        applyBestParams(*r,
                        *memory,
                        n_sources,
                        deinterlaced,
                        work,
                        max2,
                        maxVectorSizeFromBlockSizeHypothesis(max2),
                        sampleRate,
                        os,
                        args...);
        result.optimizationReport = os.str();
        {
          auto res = results.try_emplace(max2, result);
          assert(res.second);
        }
      }
      catch (std::exception const & e) {
        auto res = results.try_emplace(max2, e.what());
        assert(res.second);
        r.reset();
      }

      if(r)
      {
        auto res = reverbsByBlockSize.try_emplace(max2, std::move(r));
        assert(res.second);

        memories.push_back(std::move(memory));
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

  bool isEmpty() const {
    return reverbsByBlockSize.empty();
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
        deprecated.emplace_back(prevCurrent->reverbIt,
                                reverbUnpaddedLength);
      }
    }
    return current ? current->blockSizeHypothesys() : -1;
  }

  bool hasCurrentReverb() const {
    return static_cast<bool>(current);
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
                                                                  maxVectorSizeFromBlockSizeHypothesis(current->blockSizeHypothesys()));
      // by design, deprecated is non empty only if current is non empty
      for(auto & d : deprecated) {
        int const nApply = std::min(d.numFramesToGo, nFramesToCompute);
        d.numFramesToGo -= nApply;
        has_step_errors = !d.getReverb().addWetInputZeroVectorized(work_buffers,
                                                                   nWorkBuffers,
                                                                   nApply,
                                                                   maxVectorSizeFromBlockSizeHypothesis(d.blockSizeHypothesys())) || has_step_errors;
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

    int hostBlockSize; // may not be a power of two
    ConstIterator reverbIt;
  };
  Map reverbsByBlockSize;
  std::optional<Current> current;
  struct DeprecatedReverb {
    DeprecatedReverb(ConstIterator it, int numFramesToGo)
    : reverbIt(it),
    numFramesToGo(numFramesToGo)
    {}

    ConstIterator reverbIt;
    int numFramesToGo;

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
  };
  std::vector<DeprecatedReverb> deprecated;
  int reverbUnpaddedLength = 0;

  WorkCplxFreqs work;
  std::vector<std::unique_ptr<typename Rev::MemResource::type>> memories;

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

} // NS imajuscule::audio
