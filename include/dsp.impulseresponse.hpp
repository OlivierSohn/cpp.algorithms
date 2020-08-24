
namespace imajuscule::audio {

template<typename T>
struct DeinterlacedBuffers {
    DeinterlacedBuffers(InterlacedBuffer const & interlaced) :
    DeinterlacedBuffers(interlaced.getBuffer(), interlaced.countChannels()) {}

    template<typename Vector>
    DeinterlacedBuffers(Vector const & interlaced, int const n_channels) :
    deinterlaced(n_channels)
    {
        static_assert(std::is_same_v<typename Vector::value_type, T>);
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

  template<typename Vector>
  void deinterlace(Vector const & ir, int const n_channels) {
    static_assert(std::is_same_v<typename Vector::value_type, T>);
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
      using CplxFreqs = typename fft::RealFBins_<Tag, T, a64::Alloc>::type;

      auto N = ceil_power_of_two(v.size());

      Algo fft_algo(Contexts::getInstance().getBySize(N));
      CplxFreqs fft_of_coeffs;
      fft_of_coeffs.resize(N);
      // make a copy when passing by value
      auto coeffVec = fft::RealSignal_<Tag, T>::make(v);
      coeffVec.resize(N);
      Assert(N == coeffVec.size());
      fft_algo.forward(coeffVec.begin(), fft_of_coeffs.data(), N);
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

}
