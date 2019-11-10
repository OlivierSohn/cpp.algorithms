
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
        std::ostringstream os;
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
        std::ostringstream os;
        os << "nEarlyCofficients : " << nEarlyCofficients << std::endl;
        os << "totalSizePadded : " << totalSizePadded << std::endl;
        os << "scaleSize : " << scaleSize << std::endl;
        os << "countScales : " << countScales << std::endl;
        return os.str();
    }
};

  /*
   Depending on the number of sources, represents a convolution reverb
   or a spatialization.
   */
  template<int nAudioOut>
  struct Reverbs {
    static constexpr auto nOut = nAudioOut;
      static_assert(nAudioOut > 0);

    using ConvolutionReverb = ZeroLatencyScaledFineGrainedPartitionnedConvolution<double>;
    using Spatializer = audio::Spatializer<nAudioOut, ConvolutionReverb>;

    using FPT = typename ConvolutionReverb::FPT;
    using SetupParam = typename ConvolutionReverb::SetupParam;
    using PS = typename PartitionAlgo<ConvolutionReverb>::PS;

    int countScales() {
      if(nAudioOut && !conv_reverbs.empty() && !conv_reverbs[0].isZero()) {
        return imajuscule::countScales(conv_reverbs[0]);
      }
      else if(!spatializer.empty()) {
        return spatializer.countScales();
      }
      return 0;
    }

      /*
      // in-place, with dry/wet (see also addWet)
    void apply(double*buffer, int nFrames) {
      // raw convolutions and spatializer are mutually exclusive.
      if(nAudioOut && !conv_reverbs[0].isZero()) {
        Assert(conv_reverbs.size() == nAudioOut);
        Assert(spatializer.empty());
        for(int i=0; i<nFrames; ++i) {
          double const wet = wetRatio.step();
          double const dry = 1.-wet;
#ifndef NDEBUG
          Assert(checkDryWet(dry,wet));
#endif
          for(int j=0; j<nAudioOut; ++j) {
            auto & conv_reverb = conv_reverbs[j];
            auto & sample = buffer[i*nAudioOut + j];
            Assert(dry == 0 || conv_reverb.getLatency() == 0); // else dry and wet signals are out of sync
            sample = dry * sample + wet * conv_reverb.step(sample);
          }
        }
      }
      else if(!spatializer.empty()) {
        for(int i=0; i<nFrames; ++i) {
          double const wet = wetRatio.step();
          double const dry = 1.-wet;
#ifndef NDEBUG
          Assert(checkDryWet(dry,wet));
#endif
          spatializer.step(&buffer[i*nAudioOut], dry, wet);
        }
      }
    }
       */
      
      // out of place
      void addWet(double const * const input_buffer, double * output_buffer, int nFrames) {
          // raw convolutions and spatializer are mutually exclusive.
          if(nAudioOut && !conv_reverbs[0].isZero()) {
              Assert(conv_reverbs.size() == nAudioOut);
              Assert(spatializer.empty());
              for(int i=0; i<nFrames; ++i) {
                  for(int j=0; j<nAudioOut; ++j) {
                      int const idx = i*nAudioOut + j;
                      output_buffer[idx] += conv_reverbs[j].step(input_buffer[idx]);
                  }
              }
          }
          else if(!spatializer.empty()) {
              for(int i=0; i<nFrames; ++i) {
                  int const idx = i*nAudioOut;
                  spatializer.addWet(input_buffer + idx,
                                     output_buffer + idx);
              }
          }
      }
      
      void addWetInputZero(double * output_buffer, int nFrames) {
          // raw convolutions and spatializer are mutually exclusive.
          if(nAudioOut && !conv_reverbs[0].isZero()) {
              Assert(conv_reverbs.size() == nAudioOut);
              Assert(spatializer.empty());
              for(int i=0; i<nFrames; ++i) {
                  for(int j=0; j<nAudioOut; ++j) {
                      int const idx = i*nAudioOut + j;
                      output_buffer[idx] += conv_reverbs[j].step(0.);
                  }
              }
          }
          else if(!spatializer.empty()) {
              std::array<double, nOut> zeros{};
              for(int i=0; i<nFrames; ++i) {
                  int const idx = i*nAudioOut;
                  spatializer.addWetInputZero(output_buffer + idx);
              }
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
      for(auto & r : conv_reverbs) {
        r.reset();
      }
      spatializer.clear();
    }

    static constexpr auto ratio_hard_limit = 1.0f;
    //    because of overhead due to os managing audio, because of "other things running on the device", etc...
    // at 0.38f on ios release we have glitches when putting the app in the background
    // at 0.25f on linux we have glitches
    static constexpr auto ratio_soft_limit = 0.15f * ratio_hard_limit;

      void setConvolutionReverbIR(DeinterlacedBuffers<FPT> const & deinterlaced,
                                  int n_audiocb_frames,
                                  double sampleRate,
                                  ResponseTailSubsampling rts,
                                  std::ostream & os,
                                  ResponseStructure & structure)
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
                                                       nAudioOut);

      assert(nAudioOut == conv_reverbs.size());

      double theoretical_max_ns_per_frame = 1e9/sampleRate;

      auto mayRes =
      algo.optimize_reverb_parameters(theoretical_max_ns_per_frame * ratio_soft_limit / static_cast<float>(deinterlaced.countChannels()), rts, os);
      if(!mayRes) {
          std::stringstream ss;
          ss << "could not optimize (1)";
          throw std::runtime_error(ss.str());
      }

      auto [partitionning,n_scales] = *mayRes;
      if(!partitionning.cost) {
          std::stringstream ss;
          ss << "could not optimize (2)";
          throw std::runtime_error(ss.str());
      }
        
      setCoefficients(partitionning,
                      n_scales,
                      deinterlaced.getBuffers(),
                      structure);

      logReport(deinterlaced.countChannels(), partitionning, theoretical_max_ns_per_frame, sampleRate, os);
      algo.logReport(n_scales,partitionning, os);
    }
      
    bool hasSpatializer() const { return !spatializer.empty(); }
      
      int getUnpaddedSize() const {
          return total_response_size;
      }
      int getPaddedSize() const {
          return total_response_size_padded;
      }

  private:

    std::array<ConvolutionReverb, nAudioOut> conv_reverbs;
    Spatializer spatializer;
    int total_response_size = 0;
    int total_response_size_padded = 0;

    void setCoefficients(PS const & spec,
                         int n_scales,
                         std::vector<a64::vector<double>> deinterlaced_coeffs,
                         ResponseStructure & structure) {
      using namespace std;
      assert(n_scales >= 1);
        
        total_response_size = deinterlaced_coeffs.empty() ? 0 : deinterlaced_coeffs[0].size();
        for(auto const & v:deinterlaced_coeffs) {
            if(total_response_size != v.size()) {
                throw std::runtime_error("deinterlaced coefficients have different sizes");
            }
        }

      // debugging
      /*
       for(auto const & v : deinterlaced_coeffs) {
       StringPlot plot(40, 100);
       plot.draw(v);
       plot.log();
       }
       {
       using namespace audio;
       write_wav("/Users/Olivier/Dev/Audiofiles", "deinterlaced.wav", deinterlaced_coeffs, SAMPLE_RATE);
       }
       */



      int const n_coeffs_early_handler = std::max(minLatencyLateHandlerWhenEarlyHandlerIsDefaultOptimizedFIRFilter,
                                                  lateHandlerLatency<ConvolutionReverb>(spec.cost->partition_size));
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

      structure.totalSizePadded = total_response_size_padded;

      auto const n_channels = deinterlaced_coeffs.size();
      auto const nSources = n_channels / nAudioOut;
        
      // if we have enough sources, we can spatialize them, i.e each ear will receive
      // the sum of multiple convolutions.
      if(nSources <= 1) {
        auto n = 0;

        for(auto & rev : conv_reverbs)
        {
          prepare(*spec.cost, rev, n_scales, scale_sz);
          rev.setCoefficients(deinterlaced_coeffs[n%n_channels]);
          assert(rev.isValid());
          assert(scalesAreValid(n_scales, rev));
          assert(imajuscule::countScales(rev) == n_scales);
          // uncomment to debug
          /*
          while(!scalesAreValid(n_scales, rev)) {
            rev.reset();
            prepare(*spec.cost, rev, n_scales, scale_sz);
            rev.setCoefficients(deinterlaced_coeffs[n%n_channels]);
          }*/
          // to "spread" the computations of each channel's convolution reverbs
          // on different audio callback calls, we separate them as much as possible using a phase:
          dephase(n * spec.cost->phase, n_scales, rev);
          ++n;
        }
      }
      else {
        if(nSources * nAudioOut != n_channels) {
          throw std::logic_error("wrong number of channels");
        }

        assert(spatializer.empty());
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

          spatializer.addSourceLocation(std::move(a), *spec.cost, n_scales, scale_sz);
          assert(!spatializer.empty());
        }
        assert(!spatializer.empty());
        spatializer.dephaseComputations(spec.cost->phase, n_scales);
      }
        
    }

    template<typename T>
    int countDecimalNumbersBeforeTheDot(T val){
        val = std::abs(val);
        int count = 0;
        while(true) {
            if(val < 1) {
                return count;
            }
            val /= 10;
            ++count;
        }
    }
    void logReport(int n_channels, PS & partitionning, double theoretical_max_avg_time_per_frame, double sampleRate, std::ostream & os) {
      using namespace std;
        
        auto actual = partitionning.getCost();
        auto theoretical = theoretical_max_avg_time_per_frame / static_cast<float>(n_channels);
        auto ratio = actual / theoretical;
        
        /*
        os << "Dropouts probability    : ";
        static_assert(ratio_soft_limit < ratio_hard_limit);
        if(ratio >= ratio_soft_limit) {
            if(ratio > ratio_hard_limit) {
                os << "100 %";
            }
            else {
                os << " 50 %";
            }
        }
        else {
            os << " 0 %";
        }
        os << endl;
         */
        
        auto nActual = countDecimalNumbersBeforeTheDot(actual);
        auto nTheoretical = countDecimalNumbersBeforeTheDot(theoretical);
        std::string prefixActual =
        std::string(std::max(0, nTheoretical-nActual), ' ');
        std::string prefixTheoretical =
        std::string(std::max(0, nActual-nTheoretical), ' ');
        os << "Foreseen CPU load (1 core)          : " << std::fixed << std::setprecision(2) << 100.*ratio << "%" << endl;
        os << "Average computation time per sample : " << std::fixed << std::setprecision(0) << actual << " ns" << endl;
        //os << "- Actual                : " << prefixActual      << actual      << " ns" << endl;
        //os << "- Allowed (theoretical) : " << prefixTheoretical << theoretical << " ns" << endl;

      auto index = 1;
      for(auto const & r : conv_reverbs)
      {
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

      //os << endl;
      //os << *partitionning.cost << endl;
      constexpr auto debug_gradient_descent = false;
      if(debug_gradient_descent) {
        os << "Gradient descent report :" << endl;
        partitionning.gd.debug(true, os);
      }
    }

    void reset() {
      disable();
    }
  };

}
