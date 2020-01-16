
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

// depends on the optimization
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
        Realtime_Synchronous_Subsampled,
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
        case ReverbType::Realtime_Synchronous_Subsampled :
            return "Realtime_Synchronous_Subsampled";
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
      
    using ConvolutionReverb =
      std::conditional_t< reverbType==ReverbType::Offline,
        AlgoOptimizedFIRFilter<double, Tag>,

      std::conditional_t< reverbType==ReverbType::Realtime_Synchronous,
        AlgoZeroLatencyScaledFineGrainedPartitionnedConvolution<double, Tag>,

      std::conditional_t< reverbType==ReverbType::Realtime_Synchronous_Subsampled,
        ZeroLatencyScaledFineGrainedPartitionnedSubsampledConvolution<double, Tag>,

      std::conditional_t< reverbType==ReverbType::Realtime_Asynchronous_Legacy,
        AlgoZeroLatencyScaledAsyncConvolution<double, Tag, OnWorkerTooSlow>,

      std::conditional_t< reverbType==ReverbType::Realtime_Asynchronous,
        AlgoZeroLatencyScaledAsyncConvolutionOptimized<double, Tag, OnWorkerTooSlow>,
      void >>>>>;
      
    using Spatializer = audio::Spatializer<nAudioOut, ConvolutionReverb>;

    using FPT = typename ConvolutionReverb::FPT;
    using SetupParam = typename ConvolutionReverb::SetupParam;

    int countScales() {
        if(nAudioOut && !conv_reverbs.empty() && !conv_reverbs[0]->isZero()) {
            if constexpr(ConvolutionReverb::has_subsampling) {
                return imajuscule::countScales(*conv_reverbs[0]);
            }
            else {
                return 1;
            }
        }
        else if(!spatializer.empty()) {
            return spatializer.countScales();
        }
        return 0;
    }
      template<typename F>
      void foreachConvReverb(F f) const {
          if(nAudioOut && !conv_reverbs.empty()) {
              for(auto const & c : conv_reverbs) {
                  f(*c);
              }
          }
          else if(!spatializer.empty()) {
              return spatializer.foreachConvReverb(f);
          }
      }
      
      template<typename FPT2>
      void assignWet(FPT2 const * const * const input_buffers,
                     int nInputBuffers,
                     FPT2 ** output_buffers,
                     int nOutputBuffers,
                     int nFramesToCompute) {
          // raw convolutions and spatializer are mutually exclusive.
          if(!conv_reverbs.empty()) {
              Assert(conv_reverbs.size() == nAudioOut);
              Assert(spatializer.empty());
              Assert(nOutputBuffers == nAudioOut);
              Assert(nInputBuffers == nAudioOut);
              for(int j=0; j<nAudioOut; ++j) {
                  conv_reverbs[j]->stepAssignVectorized(input_buffers[j],
                                                        output_buffers[j],
                                                        nFramesToCompute);
              }
          }
          else if(!spatializer.empty()) {
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
          bool success = true;
          assert(vectorLength > 0);
          // raw convolutions and spatializer are mutually exclusive.
          if(!conv_reverbs.empty()) {
              Assert(conv_reverbs.size() == nAudioOut);
              Assert(spatializer.empty());
              Assert(nOutputBuffers == nAudioOut);
              Assert(nInputBuffers == nAudioOut);
              for(int j=0; j<nAudioOut; ++j) {
                  auto input_buffer = input_buffers[j];
                  auto output_buffer = output_buffers[j];
                  for(int i=0; i<nFramesToCompute; i += vectorLength) {
                    conv_reverbs[j]->stepAssignVectorized(input_buffer + i,
                                                                    output_buffer+ i,
                                                                    std::min(vectorLength, nFramesToCompute-i));
                  }
                if constexpr (ConvolutionReverb::step_can_error) {
                    success = !conv_reverbs[j]->hasStepErrors() && success;
                }
              }
          }
          else if(!spatializer.empty()) {
            success = spatializer.assignWetVectorized(input_buffers,
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
          }
        return success;
      }
      
      template<typename FPT2>
      bool addWetInputZeroVectorized(FPT2 ** output_buffers,
                                     int nOutputBuffers,
                                     int nFramesToCompute,
                                     int vectorLength) {
          bool success = true;
          assert(vectorLength > 0);
          // raw convolutions and spatializer are mutually exclusive.
          if(!conv_reverbs.empty()) {
              Assert(conv_reverbs.size() == nAudioOut);
              Assert(spatializer.empty());
              Assert(nOutputBuffers == nAudioOut);
              for(int j=0; j<nAudioOut; ++j) {
                  auto output_buffer = output_buffers[j];
                  for(int i=0; i<nFramesToCompute; i += vectorLength) {
                      conv_reverbs[j]->stepAddInputZeroVectorized(output_buffer + i,
                                                                            std::min(vectorLength, nFramesToCompute-i));
                  }
                  if constexpr (ConvolutionReverb::step_can_error) {
                    success = !conv_reverbs[j]->hasStepErrors() && success;
                  }
              }
          }
          else if(!spatializer.empty()) {
              success = spatializer.addWetInputZeroVectorized(output_buffers,
                                                              nOutputBuffers,
                                                              nFramesToCompute,
                                                              vectorLength);
          }
          return success;
      }
      
    void flushToSilence()
    {
        for(auto & r : conv_reverbs) {
            r->flushToSilence();
        }
        spatializer.flushToSilence();
    }

    void disable()
    {
      conv_reverbs.clear();
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
          
          os << "Partitioning:" << std::endl;
          
          
          auto indent = std::make_unique<IndentingOStreambuf>(os);
          using LegacyReverb = corresponding_legacy_dsp_t<ConvolutionReverb>;
          auto partitionning = PartitionAlgo<LegacyReverb>::run(n_response_channels,
                                                                nAudioOut,
                                                                n_audiocb_frames,
                                                                deinterlaced.countFrames(),
                                                                sampleRate,
                                                                max_avg_time_per_sample,
                                                                os,
                                                                args...);
          indent.reset();
          
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
          
          os << "Reports:" << std::endl;
          IndentingOStreambuf i(os);
          
          logReport(sampleRate, os);
      }
      
    bool hasSpatializer() const { return !spatializer.empty(); }
      
      // This is the 'relevant' length of the reverb: after that coefficients are all zeroes
      int getUnpaddedSize() const {
          return total_response_size;
      }
      int getPaddedSize() const {
          return total_response_size_padded;
      }

  private:

    std::vector<std::unique_ptr<ConvolutionReverb>> conv_reverbs;
    Spatializer spatializer;
    int total_response_size = 0;
    int total_response_size_padded = 0;

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

          if constexpr (ConvolutionReverb::has_subsampling) {
              int const n_scales = count_scales(spec);
              assert(n_scales >= 1);
              int lateHandlerFirstScalePartitionSize = spec.bParams.aParams.partition_size;
              int const n_coeffs_early_handler = std::max(minLatencyLateHandlerWhenEarlyHandlerIsDefaultOptimizedFIRFilter,
                                                          lateHandlerLatency<ConvolutionReverb>(lateHandlerFirstScalePartitionSize));
              int const late_response_sz = std::max(0,total_response_size - n_coeffs_early_handler);
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
        
      assert(conv_reverbs.empty());
      assert(spatializer.empty());

      // if we have more channels than sources, each ear will receive
      // the sum of multiple convolutions.
      if(n_sources == n_channels) {
          conv_reverbs.reserve(nAudioOut);
          for(int i=0; i<nAudioOut; ++i) {
              conv_reverbs.emplace_back(std::make_unique<ConvolutionReverb>());
          }

        auto n = 0;
        for(auto & prev : conv_reverbs)
        {
          auto & rev = *prev;
          rev.setup(spec);
          rev.setCoefficients(deinterlaced_coeffs[n%n_channels]);
          assert(rev.isValid());
          // uncomment to debug
          /*
          assert(scalesAreValid(n_scales, rev));
          assert(imajuscule::countScales(rev) == n_scales);
          while(!scalesAreValid(n_scales, rev)) {
            rev.reset();
            prepare(*spec, rev, n_scales, scale_sz);
            rev.setCoefficients(deinterlaced_coeffs[n%n_channels]);
          }*/
          // to "spread" the computations of each channel's convolution reverbs
          // on different audio callback calls, we separate them as much as possible using a phase:
    // TODO reenable later
          //dephase(nAudioOut, n, spec, rev);
          ++n;
        }
      }
      else {

        for(int i=0; i<n_sources; ++i) {
          std::array<a64::vector<double>, nAudioOut> a;
          for(int j=0; j<nAudioOut; ++j) {
            // for wir files of wave, it seems the order is by "ears" then by "source":

            // ear Left source 1
            // ...
            // ear Left source N
            // ear Right source 1
            // ...
            // ear Right source N

            a[j] = std::move(deinterlaced_coeffs[nAudioOut*i+j]);
          }

          spatializer.addSourceLocation(std::move(a), spec);
          assert(!spatializer.empty());
        }
        assert(!spatializer.empty());
        spatializer.dephaseComputations(spec);
      }
        
    }

    void logReport(double sampleRate, std::ostream & os) {
      using namespace std;
        
        int i=0;
        auto logConv = [&os, &i](auto const & r) mutable {
            ++i;
            os << "Convolution " << i << ":" << std::endl;
            IndentingOStreambuf indent(os);
            r.logComputeState(os);
        };

        for(auto const & pr : conv_reverbs)
        {
            logConv(*pr);
        }

        spatializer.foreachConvReverb([&logConv](auto const & r){ logConv(r); });

        if(!spatializer.empty()) {
          os <<  "Spatialization with '" << spatializer.countSources() << "' sources" << endl;
        }
        
    }

    void reset() {
      disable();
    }
  };

}
