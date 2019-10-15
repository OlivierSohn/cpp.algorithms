
namespace imajuscule
{
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
  
  static inline range<int> getScaleCountRanges(ResponseTailSubsampling rts) {
    range<int> r;
    
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
    
    return r;
  }
  
  template<typename ConvolutionReverb>
  struct ImpulseResponseOptimizer {
    using PartitionAlgo = PartitionAlgo<ConvolutionReverb>;
    using PS = typename PartitionAlgo::PS;
    
      ImpulseResponseOptimizer(audio::InterlacedBuffer const & ib, int n_audiocb_frames, int nAudioOut) :
    n_channels(ib.countChannels()),
    n_audiocb_frames(n_audiocb_frames),
    nAudioOut(nAudioOut),
    deinterlaced(ib.countChannels())
    {
      deinterlace(ib.getBuffer());
      
      truncate_irrelevant_head();
    }
    
    int countChannels() const { return n_channels; }
    
  private:
    std::vector<a64::vector<double>> deinterlaced;
    int n_channels, n_audiocb_frames, nAudioOut;
    double max_avg_time_per_sample;
    
    auto size() const { return deinterlaced.empty()? 0 : deinterlaced[0].size(); }
    
    void deinterlace(std::vector<double> const & ir) {
      auto sz = ir.size() / n_channels;
      Assert(sz * n_channels == ir.size());
      for(auto & v : deinterlaced) {
        v.reserve(sz);
      }
      for(auto it = ir.begin(), end = ir.end(); it != end;) {
        for(int i=0; i<n_channels; ++i) {
          deinterlaced[i].push_back(*it);
          ++it;
        }
      }
    }
    
    void truncate_irrelevant_head() {
      using namespace std;
      
      auto start = size();
      for(auto const & v : deinterlaced) {
        constexpr auto sliding_average_size = 15;
        constexpr auto relevant_level = .1f;
        
        auto it = find_relevant_start_relaxed(v.begin(),
                                              v.end(),
                                              relevant_level,
                                              sliding_average_size);
        auto dist = distance(v.begin(), it);
        start = min(start, static_cast<size_t>(dist));
      }
      
      if(start < size()) {
        LG(INFO, "truncating %d first sample(s)", start);
        for(auto & v : deinterlaced) {
          v.erase(v.begin(), v.begin() + start);
          v.shrink_to_fit();
        }
      }
      else {
        LG(INFO, "no truncation");
      }
    }
  public:
    
    auto &editDeinterlaced() {
      return deinterlaced;
    }
    
    Optional<std::pair<PS,int>> optimize_reverb_parameters(double max_avg_time_per_sample, ResponseTailSubsampling rts, std::ostream & os) const {
      using namespace std;
      
      auto const size_impulse_response = size();
      
      Optional<std::pair<PS,int>> res;
      
      range<int> const scales = getScaleCountRanges(rts);

      for(int n_scales = scales.getMin(); n_scales <= scales.getMax(); ++n_scales) {
        
        auto partit = PartitionAlgo::run(n_channels,
                                         n_audiocb_frames,
                                         size_impulse_response,
                                         n_scales,
                                         os);
        auto & part = partit.getWithSpread();
        if(!part.cost) {
          LG(INFO, "discarding n_scales %d", n_scales);
          continue;
        }
        
        if(!res || res->first.getCost() > part.getCost()) {
          res = {part, n_scales};
        }

        assert(res);
        if(res->first.getCost() < max_avg_time_per_sample) {
          LG(INFO, "Optimization criteria met with %d scaling levels.", n_scales);
          break;
        }
        LG(INFO, "cost %f >= %f", res->first.getCost(), max_avg_time_per_sample);
        if(n_scales == scales.getMax()) {
            throw std::runtime_error("Optimization criteria not met, there are not enough scaling levels.");
        }
      }
      
      return res;
    }
    
    void symmetrically_scale(double&scale) {
      using namespace std;
      scale = 1.;
      double avg = 0.;
      for(auto & v : deinterlaced) {
        using namespace imajuscule::fft;
        using namespace imajuscule::fft::slow_debug;
        
        using Tag = Fastest;
        
        using ScopedContext = ScopedContext_<Tag, double>;
        using Algo = Algo_<Tag, double>;
        using CplxFreqs = typename fft::RealFBins_<Tag, double>::type;
        
        auto N = ceil_power_of_two(v.size());
        
        ScopedContext setup(N);
        
        Algo fft_algo(setup.get());
        CplxFreqs fft_of_coeffs;
        fft_of_coeffs.resize(N);
        // make a copy when passing by value
        auto coeffVec = fft::RealSignal_<Tag, double>::make(v);
        coeffVec.resize(N);
        Assert(N == coeffVec.size());
        fft_algo.forward(coeffVec.begin(), fft_of_coeffs, N);
        auto unwrapped_fft_of_coeffs = unwrap_frequencies<Tag>(fft_of_coeffs, N);
        for(auto &e : unwrapped_fft_of_coeffs) {
          e *= 1 / Algo::scale;
        }
        
        auto it = unwrapped_fft_of_coeffs.begin();
        auto end = unwrapped_fft_of_coeffs.end();
        // the second half is the conjugate of the first half, so we use the first half only.
        auto mid = it + (1 + std::distance(it, end)) / 2;
        
        a64::vector<double> norms;
        norms.reserve(std::distance(it, mid));
        
        for(; it < mid; ++it) {
          // TODO avoid overrepresenting frequencies: we could group (average)
          // by frequency band (f, f*1.05)
          norms.push_back(abs(*it));
        }
        
        std::sort(norms.begin(), norms.end());
        
        Assert(!norms.empty());
        
        auto median = *(norms.begin() + norms.size() / 2);
        auto mean = std::accumulate(norms.begin(), norms.end(), 0.) / static_cast<double>(norms.size());
        //os << "avg: " << mean << " median: " << median << std::endl;
        if(mean > avg) {
          avg = mean;
        }
        
        // TODO maybe ignore too low / too high freqs ?
      }
      
      if(avg == 0) {
        return;
      }
      scale = 1. / avg;
      for(auto & v : deinterlaced) {
        for(auto &s : v) {
          s *= scale;
        }
      }
    }
    
    void logReport(int n_scales, PS const  & partitionning, std::ostream & os) const {
      using namespace std;
        os << "Render block size : " << n_audiocb_frames << " frames" << endl;
      os << "- using max fft size = " << partitionning.cost->partition_size << ", dephasing computation schedule over " << n_channels << " channel(s)." << endl;
      os << "- using ";
      // TODO range of subsampling regions: highest quality region: xxx samples / 2-subsampled region : xxx samples / 4-subsampled
        if(1 == n_scales) {
            os << "full tail resolution";
        }
        else {
            os << "reduced tail resolution with " << n_scales - 1 << " subsampling regions";
        }
        os << "." << endl;
    }
  };

}
