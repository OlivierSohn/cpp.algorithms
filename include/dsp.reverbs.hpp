
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
  
  /*
   Depending on the number of sources, represents a convolution reverb
   or a spatialization.
   */
  template<int nAudioOut>
  struct Reverbs {
    static constexpr auto nOut = nAudioOut;
    
    using ConvolutionReverb = ZeroLatencyScaledFineGrainedPartitionnedConvolution<double>;
    using Spatializer = audio::Spatializer<nAudioOut, ConvolutionReverb>;
    
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
  
    [[nodiscard]] bool setConvolutionReverbIR(std::vector<double> ir, int n_channels, int n_audiocb_frames, double sampleRate, ResponseTailSubsampling rts)
    {
      disable();
      if(n_audiocb_frames <= 0) {
        LG(WARN, "disabling reverb (%d callback size)", n_audiocb_frames);
        return false;
      }
      if(ir.size() < n_channels) {
        LG(WARN, "disabling reverb (only %d coefficients are available)", ir.size());
        return false;
      }
      using namespace std;
      
      ImpulseResponseOptimizer<ConvolutionReverb> algo(move(ir),
                                                       n_channels,
                                                       n_audiocb_frames,
                                                       nAudioOut);
      
      assert(nAudioOut == conv_reverbs.size());
      
      double theoretical_max_ns_per_frame = 1e9/sampleRate;

      auto mayRes =
      algo.optimize_reverb_parameters(theoretical_max_ns_per_frame * ratio_soft_limit / static_cast<float>(n_channels), rts);
      if(!mayRes) {
        LG(WARN, "disabling reverb (could not optimize 1)");
        return false;
      }
      
      auto [partitionning,n_scales] = *mayRes;
      if(!partitionning.cost) {
        LG(WARN, "disabling reverb (could not optimize 2)");
        return false;
      }

      algo.symmetrically_scale();
            
      setCoefficients(partitionning,
                      n_scales,
                      move(algo.editDeinterlaced()));

      algo.logReport(n_scales,partitionning);
      logReport(algo.countChannels(), partitionning, theoretical_max_ns_per_frame, sampleRate);
      return true;
    }
    
    void transitionConvolutionReverbWetRatio(double wet, int duration) {
      wetRatio.smoothly(clamp_ret(wet,0.,1.), duration);
    }
    
    bool hasSpatializer() const { return !spatializer.empty(); }
    
  private:
    
    smoothVar<double> wetRatio = {1};
    
    std::array<ConvolutionReverb, nAudioOut> conv_reverbs;
    Spatializer spatializer;
    
    void setCoefficients(PS const & spec,
                         int n_scales,
                         std::vector<a64::vector<double>> deinterlaced_coeffs) {
      using namespace std;
      assert(n_scales >= 1);
      
      // debugging
      /*
       for(auto const & v : deinterlaced_coeffs) {
       StringPlot plot(40, 100);
       plot.draw(v);
       plot.log();
       }
       {
       using namespace audio;
       write_wav("/Users/Olivier/Dev/Audiofiles", "deinterlaced.wav", deinterlaced_coeffs);
       }
       */
      
      
      
      int const n_coeffs_early_handler = std::max(minLatencyLateHandlerWhenEarlyHandlerIsDefaultOptimizedFIRFilter,
                                                  lateHandlerLatency<ConvolutionReverb>(spec.cost->partition_size));
      int const total_response_size = deinterlaced_coeffs.empty() ? 0 : deinterlaced_coeffs[0].size();
      int const late_response_sz = std::max(0,total_response_size - n_coeffs_early_handler);
      int const scale_sz = SameSizeScales::get_scale_sz(late_response_sz, n_scales);
      
      cout << "scale size: " << scale_sz << endl;
      
      if(n_scales > 1) {
        // pad the coefficients so that all scales have the same rythm.
        int const target_late_response_sz = SameSizeScales::get_max_response_sz(n_scales, scale_sz);
        for(auto & v : deinterlaced_coeffs) {
          v.resize(n_coeffs_early_handler + target_late_response_sz);
        }
      }
      
      auto const n_channels = deinterlaced_coeffs.size();
      auto nSources = n_channels / nAudioOut;
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
    
    void logReport(int n_channels, PS & partitionning, double theoretical_max_avg_time_per_frame, double sampleRate) {
      using namespace std;
      cout << "[per sample, with 0-overhead hypothesis] allowed computation time : "
      << theoretical_max_avg_time_per_frame / static_cast<float>(n_channels)
      << " ns" << endl;
      
      cout << "[avg per sample] reverb time computation : "
      << partitionning.getCost()
      << " ns" << endl;
      
      auto ratio = partitionning.getCost() * static_cast<float>(n_channels) / static_cast<float>(theoretical_max_avg_time_per_frame);
      
      cout << "ratio : " << ratio;
      static_assert(ratio_soft_limit < ratio_hard_limit);
      cout << " which will";
      if(ratio >= ratio_soft_limit) {
        if(ratio < ratio_hard_limit) {
          cout << " probably";
        }
      }
      else {
        cout << " likely not";
      }
      cout << " generate audio dropouts." << endl;
      auto index = 1;
      for(auto const & r : conv_reverbs)
      {
        if(r.isZero()) {
          continue;
        }
        cout << " reverb " << index << " :" << endl;
        {
          auto lat = r.getLatency();
          cout
          << "  latency : " << lat << " frames ("
          << lat * 1e3 / sampleRate <<  " ms)" << endl;
        }
        ++index;
      }
      
      if(!spatializer.empty()) {
        cout <<  "spatializer with '" << spatializer.countSources() << "' sources" << endl;
      }
      
      cout << *partitionning.cost << endl;
      constexpr auto debug_gradient_descent = false;
      if(debug_gradient_descent) {
        cout << "gradient descent report :" << endl;
        partitionning.gd.debug(true);
      }
    }
    
    void reset() {
      disable();
    }
  };

}
