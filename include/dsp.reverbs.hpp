
namespace imajuscule
{

#ifndef NDEBUG
inline bool checkDryWet(double dry, double wet) {
    if(wet < 0. || wet > 1.) {
        LG(ERR, "wet %f", wet);
        return false;
    }
    if(dry < 0. || dry > 1.) {
        LG(ERR, "dry %f", dry);
        return false;
    }
    return true;
}
#endif

// does not depend on the optimization
struct ResponseTopology {
    int truncatedFrames; // ignored in totalSize
    
    int totalSize; // in frames
    int nSources; // 2
    int nChannels; // 4
    
    double amplitudeNormalizationFactor;
    
    bool operator == (const ResponseTopology& o) const {
        return
        totalSize == o.totalSize &&
        nSources == o.nSources &&
        nChannels == o.nChannels;
    }
    
    std::string toString() const {
        std::stringstream os;
        os << "totalSize : " << totalSize << std::endl;
        os << "nSources : " << nSources << std::endl;
        os << "nChannels : " << nChannels << std::endl;
        return os.str();
    }
};

/*
 Depends on the optimization.
 
 It can be viewed as a "summary" of the optimized SetupParam 
 */
struct ResponseStructure {
    int nEarlyCofficients;
    int totalSizePadded;
    int scaleSize;
    int countScales;
    
    bool isConsistent() const {
        // totalSizePadded = nEarlyCofficients + scaleSize + 2*scaleSize + 4*scaleSize + 8*scaleSize
        int tmp = 0;
        int coeff = 1;
        for(int i=0; i<countScales; ++i, coeff *= 2) {
            tmp += scaleSize * coeff;
        }
        int nOverlapp = 0;
        int crossfade = scaleFadeSz::inSmallerUnits;

        for(int i=1; i<countScales; ++i, crossfade *= 2) {
            nOverlapp += crossfade;
        }
        return tmp + nEarlyCofficients - nOverlapp == totalSizePadded;
    }
    
    bool operator == (const ResponseStructure& o) const {
        return
        nEarlyCofficients==o.nEarlyCofficients &&
        totalSizePadded == o.totalSizePadded &&
        scaleSize == o.scaleSize &&
        countScales == o.countScales;
    }

    std::string toString() const {
        std::stringstream os;
        os << "nEarlyCofficients : " << nEarlyCofficients << std::endl;
        os << "totalSizePadded : " << totalSizePadded << std::endl;
        os << "scaleSize : " << scaleSize << std::endl;
        os << "countScales : " << countScales << std::endl;
        return os.str();
    }
};

struct ConvReverbOptimizationReportSuccess {
    std::string optimizationReport;
    ResponseStructure structure;
    
    bool isConsistent() const {
        return structure.isConsistent();
    }
};


static std::ostream & operator << (std::ostream &ss, const ConvReverbOptimizationReportSuccess & r) {
    ss << r.optimizationReport << std::endl;
    return ss;
}

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
    
  /*
   Depending on the number of sources, represents a convolution reverb
   or a spatialization.
   */
  template<int nAudioOut, ReverbType reverbType, PolicyOnWorkerTooSlow OnWorkerTooSlow>
  struct Reverbs {
    static constexpr auto nOut = nAudioOut;
    static_assert(nAudioOut > 0);

    using Tag = fft::Fastest;
    
    template<typename TT>
    using Allocator = AlignedAllocator<TT, Alignment::CACHE_LINE>;
      
    using ConvolutionReverb =
      std::conditional_t< reverbType==ReverbType::Offline,
        AlgoOptimizedFIRFilter<double, Allocator, Tag>,

      std::conditional_t< reverbType==ReverbType::Realtime_Synchronous,
        AlgoZeroLatencyScaledFineGrainedPartitionnedConvolution<double, Allocator, Tag>,

      //std::conditional_t< reverbType==ReverbType::Realtime_Synchronous_Subsampled,
      //  ZeroLatencyScaledFineGrainedPartitionnedSubsampledConvolution<double, Allocator, Tag>,

      std::conditional_t< reverbType==ReverbType::Realtime_Asynchronous_Legacy,
        AlgoZeroLatencyScaledAsyncConvolution<double, Allocator, Tag, OnWorkerTooSlow>,

      std::conditional_t< reverbType==ReverbType::Realtime_Asynchronous,
        AlgoZeroLatencyScaledAsyncConvolutionOptimized<double, Allocator, Tag, OnWorkerTooSlow>,
    
      void
    
    >>>>/*>*/;
      
    using Spatializer = audio::Spatializer<nAudioOut, ConvolutionReverb>;

    using FPT = typename ConvolutionReverb::FPT;
    using SetupParam = typename ConvolutionReverb::SetupParam;

    int countScales() {
        return spatializer.countScales();
    }
      template<typename F>
      void foreachConvReverb(F f) const {
          spatializer.foreachConvReverb(f);
      }
      
      template<typename FPT2>
      void assignWet(FPT2 const * const * const input_buffers,
                     int nInputBuffers,
                     FPT2 ** output_buffers,
                     int nOutputBuffers,
                     int nFramesToCompute) {
          if(!spatializer.empty()) {
              spatializer.assignWet(input_buffers,
                                    nInputBuffers,
                                    output_buffers,
                                    nOutputBuffers,
                                    nFramesToCompute);
          }
          else {
              // zero output_buffer
              for(int j=0; j<nAudioOut; ++j) {
                  auto output_buffer = output_buffers[j];
                  fft::RealSignal_<fft::Fastest, FPT2>::zero_n_raw(output_buffer, nFramesToCompute);
              }
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
          if(!spatializer.empty()) {
            return spatializer.assignWetVectorized(input_buffers,
                                                       nInputBuffers,
                                                       output_buffers,
                                                       nOutputBuffers,
                                                       nFramesToCompute,
                                                       vectorLength);
          }
          else {
              // zero output_buffer
              for(int j=0; j<nAudioOut; ++j) {
                  auto output_buffer = output_buffers[j];
                  fft::RealSignal_<fft::Fastest, FPT2>::zero_n_raw(output_buffer, nFramesToCompute);
              }
              return true;
          }
      }
      
      template<typename FPT2>
      bool addWetInputZeroVectorized(FPT2 ** output_buffers,
                                     int nOutputBuffers,
                                     int nFramesToCompute,
                                     int vectorLength) {
          assert(vectorLength > 0);
          if(!spatializer.empty()) {
              return spatializer.addWetInputZeroVectorized(output_buffers,
                                                              nOutputBuffers,
                                                              nFramesToCompute,
                                                              vectorLength);
          }
          return true;
      }
      
    void flushToSilence()
    {
        spatializer.flushToSilence();
    }

    void disable()
    {
      spatializer.clear();
    }

    template<typename ...Args>
      void setConvolutionReverbIR(int const n_sources,
                                  DeinterlacedBuffers<FPT> const & deinterlaced,
                                  int n_audiocb_frames,
                                  double sampleRate,
                                  std::ostream & os,
                                  ResponseStructure & structure,
                                  Args... args)
      {
          disable();
          
          if(n_audiocb_frames <= 0) {
              std::stringstream ss;
              ss << "Negative or zero callback size (" << n_audiocb_frames << ")";
              throw std::runtime_error(ss.str());
          }
          
          using namespace std;
          
          int const n_response_channels = deinterlaced.countChannels();
          
          static constexpr double ratio_hard_limit = 1.0;
          //    because of overhead due to os managing audio, because of "other things running on the device", etc...
          // at 0.38f on ios release we have glitches when putting the app in the background
          // at 0.25f on linux we have glitches
          static constexpr double ratio_soft_limit = 0.15 * ratio_hard_limit;
          
          double const theoretical_max_ns_per_frame(1e9/sampleRate);
          double const max_avg_time_per_sample(theoretical_max_ns_per_frame * ratio_soft_limit / static_cast<float>(n_response_channels));

          using LegacyReverb = corresponding_legacy_dsp_t<ConvolutionReverb>;
          auto partitionning = PartitionAlgo<LegacyReverb>::run(n_response_channels,
                                                                nAudioOut,
                                                                n_audiocb_frames,
                                                                deinterlaced.countFrames(),
                                                                sampleRate,
                                                                max_avg_time_per_sample,
                                                                os,
                                                                args...);
          if(!partitionning) {
              std::stringstream ss;
              ss << "could not optimize (2) :" << std::endl << os.rdbuf();
              throw std::runtime_error(ss.str());
          }
          // buffers are copied here:
          auto buffers = deinterlaced.getBuffers();
          structure = handleScales(*partitionning,
                                   buffers);
          
          partitionning->logReport(deinterlaced.countChannels(),
                                   theoretical_max_ns_per_frame,
                                   os);

          setCoefficients(n_sources,
                          *partitionning,
                          buffers);
          logReport(sampleRate, os);
      }
      
      // This is the 'relevant' length of the reverb: after that coefficients are all zeroes
      int getUnpaddedSize() const {
          return total_response_size;
      }

  private:

    Spatializer spatializer;
    int total_response_size = 0;

      ResponseStructure handleScales(SetupParam const & spec,
                                     std::vector<a64::vector<double>> & deinterlaced_coeffs) {
        using namespace std;
        ResponseStructure structure;
          
          total_response_size = deinterlaced_coeffs.empty() ? 0 : deinterlaced_coeffs[0].size();
          for(auto const & v:deinterlaced_coeffs) {
              if(total_response_size != v.size()) {
                  throw std::runtime_error("deinterlaced coefficients have different sizes");
              }
          }

          int total_response_size_padded = 0;

          if constexpr (ConvolutionReverb::has_subsampling) {
              int const n_scales = count_scales(spec);
              assert(n_scales >= 1);
              int lateHandlerFirstScalePartitionSize = spec.b.a.partition_size;
              int const n_coeffs_early_handler = std::max(minLatencyLateHandlerWhenEarlyHandlerIsDefaultOptimizedFIRFilter,
                                                          earliestDeepestLatency<typename ConvolutionReverb::LateHandler>(lateHandlerFirstScalePartitionSize)).toInteger();
              int const late_response_sz = std::max(0
                                                    ,total_response_size - n_coeffs_early_handler);
              int const scale_sz = SameSizeScales::get_scale_sz(late_response_sz, n_scales);

                structure.scaleSize = scale_sz;
                structure.countScales = n_scales;
                structure.nEarlyCofficients = n_coeffs_early_handler;

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
              structure.scaleSize = 0;
              structure.countScales = 1;
              total_response_size_padded = total_response_size;
              structure.nEarlyCofficients = total_response_size_padded;
          }
          structure.totalSizePadded = total_response_size_padded;
          return structure;
      }
      
    void setCoefficients(int const n_sources,
                         SetupParam const & spec,
                         std::vector<a64::vector<double>> const & deinterlaced_coeffs) {
      using namespace std;
        
      auto const n_channels = deinterlaced_coeffs.size();
      Assert(nAudioOut <= n_channels);
        int const n_conv_per_source = n_channels / n_sources;
        if(n_channels != n_conv_per_source * n_sources) {
            throw std::runtime_error("inconsistent number of audio sources / channels");
        }
        
        assert(spatializer.empty());
        // for wir files of wave, it seems the order is by "ears" then by "source":

        // ear Left source 1
        // ...
        // ear Left source N
        // ear Right source 1
        // ...
        // ear Right source N
        spatializer.setSources(n_sources,
                               std::move(deinterlaced_coeffs),
                               spec);

        assert(!spatializer.empty());
        spatializer.dephaseComputations();
    }

    void logReport(double sampleRate, std::ostream & os)
    {
        spatializer.logReport(sampleRate, os);
    }

    void reset() {
      disable();
    }
  };

}
