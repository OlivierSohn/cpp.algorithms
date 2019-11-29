
namespace imajuscule
{

template<typename T>
struct DeinterlacedBuffers {
    
    DeinterlacedBuffers(std::vector<T> const & interlaced, int const n_channels) :
    deinterlaced(n_channels)
    {
        deinterlace(interlaced, n_channels);
    }
    
    DeinterlacedBuffers(std::vector<a64::vector<T>> deinterlaced) :
    deinterlaced(std::move(deinterlaced))
    {}
    
    int countChannels() const {
        return deinterlaced.size();
    }
    
    int countFrames() const {
        return deinterlaced.empty()? 0 : deinterlaced[0].size();
    }
    
    auto const & getBuffers() const {
        return deinterlaced;
    }
    
private:
  std::vector<a64::vector<T>> deinterlaced;
  
  void deinterlace(std::vector<T> const & ir, int const n_channels) {
    if(ir.size() < n_channels) {
        std::stringstream ss;
        ss << "Not enough coefficients : " << ir.size() << " for " << n_channels << " channel(s).";
        throw std::runtime_error(ss.str());
    }
      if(n_channels == 0) {
          return;
      }
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
public:

    // returns the size of the truncation
    int truncate_irrelevant_head() {
        using namespace std;
        
        int count = countFrames();
        for(auto const & v : deinterlaced) {
            constexpr auto sliding_average_size = 15;
            constexpr auto relevant_level = .1f;
            
            auto it = find_relevant_start_relaxed(v.begin(),
                                                  v.end(),
                                                  relevant_level,
                                                  sliding_average_size);
            auto dist = distance(v.begin(), it);
            count = min(count, static_cast<int>(dist));
        }
        
        if(count >= countFrames()) {
            // no truncation
            return 0;
        }
        // truncate first 'count' samples
        for(auto & v : deinterlaced) {
            v.erase(v.begin(), v.begin() + count);
        }
        return count;
    }
  
  auto &editDeinterlaced() {
    return deinterlaced;
  }
  
  T symmetrically_scale() {
    using namespace std;
    T avg = 0.;
    for(auto & v : deinterlaced) {
      using namespace imajuscule::fft;
      using namespace imajuscule::fft::slow_debug;
      
      using Tag = Fastest;
      
      using Contexts = fft::Contexts_<Tag, T>;
      using Algo = Algo_<Tag, T>;
      using CplxFreqs = typename fft::RealFBins_<Tag, T>::type;
      
      auto N = ceil_power_of_two(v.size());
      
      Algo fft_algo(Contexts::getInstance().getBySize(N));
      CplxFreqs fft_of_coeffs;
      fft_of_coeffs.resize(N);
      // make a copy when passing by value
      auto coeffVec = fft::RealSignal_<Tag, T>::make(v);
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
      
      a64::vector<T> norms;
      norms.reserve(std::distance(it, mid));
      
      for(; it < mid; ++it) {
        // TODO avoid overrepresenting frequencies: we could group (average)
        // by frequency band (f, f*1.05)
        norms.push_back(abs(*it));
      }
      
      std::sort(norms.begin(), norms.end());
      
      Assert(!norms.empty());
      
      auto median = *(norms.begin() + norms.size() / 2);
      auto mean = std::accumulate(norms.begin(), norms.end(), 0.) / static_cast<T>(norms.size());
      //os << "avg: " << mean << " median: " << median << std::endl;
      if(mean > avg) {
        avg = mean;
      }
      
      // TODO maybe ignore too low / too high freqs ?
    }
    
    if(avg == 0) {
      return 1.;
    }
    T scale = 1. / avg;
    for(auto & v : deinterlaced) {
      for(auto &s : v) {
        s *= scale;
      }
    }
      return scale;
  }
  
};

  template<typename ConvolutionReverb>
  struct ImpulseResponseOptimizer {
    using PartitionAlgo = PartitionAlgo<ConvolutionReverb>;
    using PS = typename PartitionAlgo::PS;
    using FPT = typename ConvolutionReverb::FPT;
      
      static constexpr auto ratio_hard_limit = 1.0f;
      //    because of overhead due to os managing audio, because of "other things running on the device", etc...
      // at 0.38f on ios release we have glitches when putting the app in the background
      // at 0.25f on linux we have glitches
      static constexpr auto ratio_soft_limit = 0.15f * ratio_hard_limit;

      ImpulseResponseOptimizer(int const n_response_channels,
                               int const n_response_frames,
                               int n_audiocb_frames,
                               int nAudioOut,
                               double frame_rate) :
      n_audiocb_frames(n_audiocb_frames),
      nAudioOut(nAudioOut),
      n_response_channels(n_response_channels),
      n_response_frames(n_response_frames),
      frame_rate(frame_rate),
      theoretical_max_ns_per_frame(1e9/frame_rate),
      max_avg_time_per_sample(theoretical_max_ns_per_frame * ratio_soft_limit / static_cast<float>(n_response_channels))
      {}
      
  private:
      int n_audiocb_frames, nAudioOut, n_response_channels, n_response_frames;
      double frame_rate;
  public:
      double const theoretical_max_ns_per_frame;
  private:
      double const max_avg_time_per_sample;
      
  public:
      
      template<typename ...Args>
      auto optimize_reverb_parameters(std::ostream & os, Args ...args) const
      {
          return PartitionAlgo::run(n_response_channels,
                                    nAudioOut,
                                    n_audiocb_frames,
                                    n_response_frames,
                                    frame_rate,
                                    max_avg_time_per_sample,
                                    os,
                                    args...);
      }
    
    void logReport(std::ostream & os) const {
      using namespace std;
        os << "Render block size : " << n_audiocb_frames << " frames" << endl;
        os << "- dephasing computation schedule over " << n_response_channels << " channel(s)" << endl;
    }
  };

}
