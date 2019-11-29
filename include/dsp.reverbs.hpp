
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

    enum class ReverbType {
        Offline,
        Realtime_Synchronous,
        Realtime_Synchronous_Subsampled,
        Realtime_Asynchronous
    };

static inline std::string toString(ReverbType t) {
    switch(t) {
        case ReverbType::Offline :
            return "Offline";
        case ReverbType::Realtime_Synchronous :
            return "Realtime_Synchronous";
        case ReverbType::Realtime_Synchronous_Subsampled :
            return "Realtime_Synchronous_Subsampled";
        case ReverbType::Realtime_Asynchronous :
            return "Realtime_Asynchronous";
    }
    return "?";
}
  /*
   Depending on the number of sources, represents a convolution reverb
   or a spatialization.
   */
  template<int nAudioOut, ReverbType reverbType>
  struct Reverbs {
    static constexpr auto nOut = nAudioOut;
      static_assert(nAudioOut > 0);

    using ConvolutionReverb =
      std::conditional_t< reverbType==ReverbType::Offline,
        OptimizedFIRFilter<double>,

      std::conditional_t< reverbType==ReverbType::Realtime_Synchronous,
        ZeroLatencyScaledFineGrainedPartitionnedConvolution<double>,

      std::conditional_t< reverbType==ReverbType::Realtime_Synchronous_Subsampled,
        ZeroLatencyScaledFineGrainedPartitionnedSubsampledConvolution<double>,

      std::conditional_t< reverbType==ReverbType::Realtime_Asynchronous,
        ZeroLatencyScaledAsyncConvolution<double>,
      void >>>>;
      
    using Spatializer = audio::Spatializer<nAudioOut, ConvolutionReverb>;

    using FPT = typename ConvolutionReverb::FPT;
    using SetupParam = typename ConvolutionReverb::SetupParam;
    using PS = typename PartitionAlgo<ConvolutionReverb>::PS;

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
                  using FFTTag = fft::Fastest;
                  fft::RealSignal_<FFTTag, FPT2>::zero_n_raw(output_buffer, nFramesToCompute);
              }
          }
      }
      template<typename FPT2>
      void assignWetVectorized(FPT2 const * const * const input_buffers,
                               int nInputBuffers,
                               FPT2 ** output_buffers,
                               int nOutputBuffers,
                               int nFramesToCompute,
                               int vectorLength) {
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
              }
          }
          else if(!spatializer.empty()) {
              spatializer.assignWetVectorized(input_buffers,
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
                  using FFTTag = fft::Fastest;
                  fft::RealSignal_<FFTTag, FPT2>::zero_n_raw(output_buffer, nFramesToCompute);
              }
          }
      }
      
      template<typename FPT2>
      void addWetInputZeroVectorized(FPT2 ** output_buffers,
                                     int nOutputBuffers,
                                     int nFramesToCompute,
                                     int vectorLength) {
          assert(vectorLength > 0);
          // raw convolutions and spatializer are mutually exclusive.
          if(nAudioOut && !conv_reverbs[0].isZero()) {
              Assert(conv_reverbs.size() == nAudioOut);
              Assert(spatializer.empty());
              Assert(nOutputBuffers == nAudioOut);
              for(int j=0; j<nAudioOut; ++j) {
                  auto output_buffer = output_buffers[j];
                  for(int i=0; i<nFramesToCompute; i += vectorLength) {
                      conv_reverbs[j].stepAddInputZeroVectorized(output_buffer + i,
                                                                 std::min(vectorLength, nFramesToCompute-i));
                  }
              }
          }
          else if(!spatializer.empty()) {
              spatializer.addWetInputZeroVectorized(output_buffers,
                                                    nOutputBuffers,
                                                    nFramesToCompute,
                                                    vectorLength);
          }
      }
      
    void flushToSilence()
    {
        for(auto & r : conv_reverbs) {
            r.flushToSilence();
        }
        spatializer.flushToSilence();
    }

    void disable()
    {
      conv_reverbs.clear();
      spatializer.clear();
    }

    template<typename ...Args>
      void setConvolutionReverbIR(DeinterlacedBuffers<FPT> const & deinterlaced,
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
      ImpulseResponseOptimizer<ConvolutionReverb> algo(deinterlaced.countChannels(),
                                                       deinterlaced.countFrames(),
                                                       n_audiocb_frames,
                                                       nAudioOut,
                                                       sampleRate);

      auto partitionning = algo.optimize_reverb_parameters(os, args...).getWithSpread();
        /*
      if(!mayRes) {
          std::stringstream ss;
          ss << "could not optimize (1) :" << std::endl << os.rdbuf();
          throw std::runtime_error(ss.str());
      }

      auto [partitionning,n_scales] = *mayRes;
      */
        if(!partitionning.cost) {
            std::stringstream ss;
            ss << "could not optimize (2) :" << std::endl << os.rdbuf();
            throw std::runtime_error(ss.str());
        }

        auto buffers = deinterlaced.getBuffers();
        auto const & param = handleScales(partitionning.cost.value(),
                        buffers,
                        structure);

      setCoefficients(param,
                      buffers);

        logReport(sampleRate, os);
        algo.logReport(os);

        partitionning.logReport(deinterlaced.countChannels(),
                                algo.theoretical_max_ns_per_frame,
                                os);
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

      template<typename U>
      SetupParam const & handleScales(U const & spec,
                           std::vector<a64::vector<double>> & deinterlaced_coeffs,
                           ResponseStructure & structure) {
        using namespace std;
          
          
          total_response_size = deinterlaced_coeffs.empty() ? 0 : deinterlaced_coeffs[0].size();
          for(auto const & v:deinterlaced_coeffs) {
              if(total_response_size != v.size()) {
                  throw std::runtime_error("deinterlaced coefficients have different sizes");
              }
          }

          SetupParam const * p = nullptr;
          if constexpr (ConvolutionReverb::has_subsampling) {
              assert(spec.n_scales >= 1);
              int lateHandlerFirstScalePartitionSize = spec.val.bParams.aParams.partition_size;
              int const n_coeffs_early_handler = std::max(minLatencyLateHandlerWhenEarlyHandlerIsDefaultOptimizedFIRFilter,
                                                          lateHandlerLatency<ConvolutionReverb>(lateHandlerFirstScalePartitionSize));
              int const late_response_sz = std::max(0,total_response_size - n_coeffs_early_handler);
              int const scale_sz = SameSizeScales::get_scale_sz(late_response_sz, spec.n_scales);

                structure.scaleSize = scale_sz;
                structure.countScales = spec.n_scales;
                structure.nEarlyCofficients = n_coeffs_early_handler;

                if(spec.n_scales > 1) {
                    // pad the coefficients so that all scales have the same rythm.
                    int const target_late_response_sz = SameSizeScales::get_max_response_sz(spec.n_scales, scale_sz);
                    
                    total_response_size_padded = n_coeffs_early_handler + target_late_response_sz;
                }
                else {
                    total_response_size_padded = total_response_size;
                }
                
                for(auto & v : deinterlaced_coeffs) {
                    v.resize(total_response_size_padded);
                }
              
              static_assert(std::is_same_v<U, WithNScales<SetupParam>>);
              p = &spec.val;
          }
          else {
              structure.scaleSize = 0;
              structure.countScales = 1;
              total_response_size_padded = total_response_size;
              structure.nEarlyCofficients = total_response_size_padded;
              
              static_assert(std::is_same_v<U, SetupParam>);
              p = &spec;
          }
          structure.totalSizePadded = total_response_size_padded;
          
          return *p;
      }
      
    void setCoefficients(SetupParam const & spec,
                         std::vector<a64::vector<double>> const & deinterlaced_coeffs) {
      using namespace std;
        
      auto const n_channels = deinterlaced_coeffs.size();
      auto const nSources = n_channels / nAudioOut;
        
      assert(conv_reverbs.empty());
      assert(spatializer.empty());

        // if we have enough sources, we can spatialize them, i.e each ear will receive
      // the sum of multiple convolutions.
      if(nSources <= 1) {
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
            prepare(*spec.cost, rev, n_scales, scale_sz);
            rev.setCoefficients(deinterlaced_coeffs[n%n_channels]);
          }*/
          // to "spread" the computations of each channel's convolution reverbs
          // on different audio callback calls, we separate them as much as possible using a phase:
          dephase(nAudioOut, n, spec, rev);
          ++n;
        }
      }
      else {
        if(nSources * nAudioOut != n_channels) {
          throw std::logic_error("wrong number of channels");
        }

        for(int i=0; i<nSources; ++i) {
          std::array<a64::vector<double>, nAudioOut> a;
          for(int j=0; j<nAudioOut; ++j) {
            // for wir files of wave, it seems the order is by "ears" then by "source":

            // ear Left source 1
            // ...
            // ear Left source N
            // ear Right source 1
            // ...
            // ear Right source N

            a[j] = std::move(deinterlaced_coeffs[i+nAudioOut*j]);
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

        auto index = 1;
        for(auto const & pr : conv_reverbs)
        {
            auto & r = *pr;
          if(r.isZero()) {
            continue;
          }
            auto lat = r.getLatency();
            if(lat) {
                os << "Channel " << index << " latency : " << lat << " frames (" << lat * 1e3 / sampleRate <<  " ms)" << endl;
            }
          ++index;
        }

        if(!spatializer.empty()) {
          os <<  "Spatialization with '" << spatializer.countSources() << "' sources" << endl;
        }
        
    }

    void reset() {
      disable();
    }
  };

}
