
namespace imajuscule
{
  template<typename ConvolutionReverb>
  struct ImpulseResponseOptimizer {
    using PartitionAlgo = PartitionAlgo<ConvolutionReverb>;
    using PS = typename PartitionAlgo::PS;
    
    ImpulseResponseOptimizer(std::vector<double> ir, int n_channels, int n_audiocb_frames, int nAudioOut) :
    n_channels(n_channels),
    n_audiocb_frames(n_audiocb_frames),
    nAudioOut(nAudioOut),
    deinterlaced(n_channels)
    {
      deinterlace(std::move(ir));
      
      truncate_irrelevant_head();
    }
    
    int countChannels() const { return n_channels; }
    
  private:
    std::vector<a64::vector<double>> deinterlaced;
    int n_channels, n_audiocb_frames, nAudioOut;
    double max_avg_time_per_sample;
    
    auto size() const { return deinterlaced.empty()? 0 : deinterlaced[0].size(); }
    
    void deinterlace(std::vector<double> ir) {
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
    
    Optional<std::pair<PS,int>> optimize_reverb_parameters(double max_avg_time_per_sample) const {
      using namespace std;
      
      auto const size_impulse_response = size();
      
      Optional<std::pair<PS,int>> res;
      
      for(int n_scales = 1; n_scales <= 4; ++n_scales) {
        
        auto partit = PartitionAlgo::run(n_channels,
                                         n_audiocb_frames,
                                         size_impulse_response,
                                         n_scales);
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
        if(n_scales == 4) {
          LG(WARN, "Optimization criteria not met, there are not enough scaling levels.");
        }
      }
      
      return res;
    }
    
    void symmetrically_scale() {
      using namespace std;
      double avg = 0.;
      for(auto & v : deinterlaced) {
        using namespace imajuscule::fft;
        using namespace imajuscule::fft::slow_debug;
        
        using Tag = Fastest;
        
        using ScopedContext = ScopedContext_<Tag, double>;
        using Algo = Algo_<Tag, double>;
        using RealSignal = typename fft::RealSignal_<Tag, double>::type;
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
        LG(INFO,"avg: %f median: %f", mean, median);
        if(mean > avg) {
          avg = mean;
        }
        
        // TODO maybe ignore too low / too high freqs ?
      }
      
      if(avg == 0) {
        return;
      }
      auto scale = 1. / avg;
      LG(INFO, "impulse response volume scaled by : %f", scale);
      for(auto & v : deinterlaced) {
        for(auto &s : v) {
          s *= scale;
        }
      }
    }
    
    void logReport(int n_scales, PS const  & partitionning) const {
      using namespace std;
      cout << endl;
      cout << "finished optimizing partition size for " << n_audiocb_frames << " cb frames" << endl;
      cout << "using " << n_scales << " scale(s)" << endl;
      cout << "using partition size " << partitionning.cost->partition_size << " with spread on " << n_channels << " channel(s)." << endl;
    }
    
  };

}
