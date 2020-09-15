namespace imajuscule::audio {

/** Returns the frequencies of a windowed signal.
 *
 * For better efficiency, the signal length should be a power of 2
 * @param input_stride stride for input
 * @param half_window is the right window half, starting at 1.0 (hence the full window has 2 times the 1.0 value)
 * @param frequencies_sqmag is the sqared magnitude of the frequencies from 0 to nyquist freq included.
 */
template<typename ITER>
void findFrequenciesSqMag(ITER it,
                     ITER const end,
                     int const windowed_signal_stride, // use to ignore high frequencies
                     std::vector<typename ITER::value_type> const & half_window,
                     int const zero_padding_factor,
                     std::vector<float> & frequencies_sqmag) {
  using namespace fft;
  using VAL = typename ITER::value_type;
  using Tag = Fastest;
  
  using Algo = Algo_<Tag, VAL>;
  using RealFBins = RealFBins_<Tag, VAL, a64::Alloc>;
  using ScopedContext = ScopedContext_<Tag, VAL>;
  
  a64::vector<VAL> v;
  int const numSamplesUnstrided = std::distance(it, end);
  int const numSamplesStrided = 1 + (numSamplesUnstrided-1) / windowed_signal_stride;
  Assert(numSamplesStrided == 2 * half_window.size());
  auto fft_length = ceil_power_of_two(numSamplesStrided * zero_padding_factor);
  v.reserve(fft_length);

  // apply the window...
  auto wEnd = half_window.end();
  auto wBegin = half_window.begin();
  auto wIt = wEnd - 1;
  bool first_half = true;
  for(; it < end; it += windowed_signal_stride) {

    Assert(wIt < wEnd);
    Assert(wIt >= wBegin);

    v.push_back(*it * *wIt);

    if (first_half) {
      if (wIt == wBegin) {
        first_half = false;
      } else {
        --wIt;
      }
    } else {
      ++wIt;
    }
  }

  // ... and pad
  v.resize(fft_length, VAL{0});

  auto signal = RealSignal_<Tag, VAL>::make(std::move(v));
  
  typename RealFBins::type result(fft_length);
  
  ScopedContext scoped_context(fft_length);
  Algo fft(scoped_context.get());
  fft.forward(signal.begin(), result.data(), fft_length);
  unwrap_frequencies_sqmag<Tag>(result, fft_length, frequencies_sqmag);

  constexpr VAL inv_scale_squared = 1. / (Algo::scale * Algo::scale);
  for (auto & v : frequencies_sqmag) {
    v *= inv_scale_squared;
  }
}

struct FreqSqMag {
  float freq; // unit: bins
  float sq_mag;
  
  bool operator == (FreqSqMag const & o) const {
    return std::make_tuple(freq, sq_mag) ==
      std::make_tuple(o.freq, o.sq_mag);
  }
};

template<typename T>
struct QuadraticInterpolation {
  QuadraticInterpolation(T const alpha, T const beta, T const gamma) {
    Assert(beta >= alpha && beta > gamma); // enforced by extractMaxFreqSqMag

    T denom = alpha - 2.*beta + gamma;
    Assert(denom >= 0);
    if (denom == 0) {
      p = 0.;
    } else {
      p = 0.5 * (alpha - gamma) / denom;
    }
    mag = beta - 0.25 * (alpha - gamma) * p;
  }
  
  auto getMag() const { return mag; }
  auto getP() const { return p; }

 private:
  T p, mag;
};

// uses quadratic interpolation
template<typename T, typename TransformAmplitude>
FreqSqMag extractMaxFreqSqMag(std::vector<T> const & freqs_sqmag, TransformAmplitude transform_amplitude) {
  Assert(!freqs_sqmag.empty());
  Optional<FreqSqMag> res;
  for (int i=0, sz=freqs_sqmag.size(); i<sz; ++i) {
    if (i < sz-1) {
      if (freqs_sqmag[i+1] >= freqs_sqmag[i]) {
        continue;
      }
      // beta > gamma (cf QuadraticInterpolation)
    }
    if (i > 0) {
      if (freqs_sqmag[i-1] > freqs_sqmag[i]) {
        continue;
      }
      // beta >= alpha (cf QuadraticInterpolation)
    }
    // i is a local maximum
    float sq_mag;
    float frequency;
    if (unlikely(i == 0 || i == sz-1)) {
      // we cannot interpolate, we are at a boundary
      sq_mag = transform_amplitude(freqs_sqmag[i]);
      frequency = i;
    } else {
      // use quadratic interpolation
      auto itp = QuadraticInterpolation(transform_amplitude(freqs_sqmag[i-1]),
                                        transform_amplitude(freqs_sqmag[i]),
                                        transform_amplitude(freqs_sqmag[i+1]));
      frequency = i + itp.getP();
      sq_mag = itp.getMag();
    }
    if (!res || res->sq_mag < sq_mag) {
      res = {frequency, sq_mag};
    }
  }
  Assert(res);
  return *res;
}

template<typename Iter>
void findMaxMagFrequencies(Iter it,
                           Iter const end,
                           int const windowed_signal_stride, // use to ignore high frequencies
                           std::vector<typename Iter::value_type> const & half_window,
                           int window_center_stride,
                           std::vector<FreqSqMag> & freqs_sqmags) {
  using T = typename Iter::value_type;
  
  int const szWindow = 2 * half_window.size();
  int const windowSampleSpan = (szWindow - 1) * windowed_signal_stride + 1;
  int const numSamples = std::distance(it, end);
  int const num_ffts = 1 + (numSamples - windowSampleSpan) / window_center_stride;
  
  auto limit = it + windowSampleSpan;
  
  freqs_sqmags.clear();
  freqs_sqmags.reserve(num_ffts);
  
  std::vector<float> freqs_sqmag;
  for (; limit <= end; it += window_center_stride, limit += window_center_stride) {
    findFrequenciesSqMag(it, limit, windowed_signal_stride, half_window, 1, freqs_sqmag);
    freqs_sqmags.push_back(extractMaxFreqSqMag(freqs_sqmag, [](T const val){
      Assert(val >= 0);
      auto const s = std::sqrt(val);
      if (s == 0.) {
        return std::numeric_limits<T>::min();
      }
      Assert(s > 0);
      return std::log(s);
    }));
  }
  Assert(freqs_sqmags.size() == num_ffts);
}

template<typename Iter>
void findAllMagFrequencies(Iter it,
                           Iter const end,
                           int const windowed_signal_stride, // use to ignore high frequencies
                           std::vector<typename Iter::value_type> const & half_window,
                           int window_center_stride,
                           int const zero_padding_factor,
                           std::vector<std::vector<float>> & freqs_mags) {
  int const szWindow = 2 * half_window.size();
  int const windowSampleSpan = (szWindow - 1) * windowed_signal_stride + 1;
  int const numSamples = std::distance(it, end);
  int const num_ffts = 1 + (numSamples - windowSampleSpan) / window_center_stride;
  
  auto limit = it + windowSampleSpan;

  freqs_mags.clear();
  freqs_mags.resize(num_ffts);
  
  auto itRes = freqs_mags.begin();
  for (; limit <= end; it += window_center_stride, limit += window_center_stride, ++itRes) {
    Assert(itRes < freqs_mags.end());
    findFrequenciesSqMag(it, limit, windowed_signal_stride, half_window, zero_padding_factor, *itRes);
  }
  Assert(freqs_mags.size() == num_ffts);
}

template<typename T>
void drawSpectrum(WAVReader & reader,
                  int const windowed_signal_stride, // use to ignore high frequencies
                  std::vector<T> const & half_window,
                  int window_center_stride,
                  int const zero_padding_factor,
                  std::string const& imageFile) {
  std::vector<std::vector<T>> deinterlaced;
  read_wav_as_floats(reader, deinterlaced);

  std::vector<std::vector<T>> freqs_mags;
  
  for (auto const & v : deinterlaced) {
    findAllMagFrequencies(v.begin(),
                          v.end(),
                          windowed_signal_stride,
                          half_window,
                          window_center_stride,
                          zero_padding_factor,
                          freqs_mags);
    break; // use only first channel
  }
  std::vector<T> all_freqs_mags;
  all_freqs_mags.reserve(freqs_mags.size() * freqs_mags.begin()->size());
  for (int i=0, sz = freqs_mags.begin()->size(); i<sz; ++i) {
    for (auto const & v : freqs_mags) {
      all_freqs_mags.push_back(v[i]);
    }
  }
  
  Optional<T> max_mag;
  for (auto const & val : all_freqs_mags) {
    Assert(val >= 0.f);
    if (!max_mag || val > *max_mag) {
      max_mag = val;
    }
  }
  
  Assert(max_mag);
  std::vector<unsigned char> data;
  data.reserve(3 * all_freqs_mags.size());
  for (auto const & val : all_freqs_mags) {
    T ratio = val / *max_mag;
    Assert(ratio >= 0.f);
    Assert(ratio <= 1.f);
    data.push_back(static_cast<unsigned char>(ratio * 255.9)); // r
    data.push_back(static_cast<unsigned char>(ratio * 255.9)); // g
    data.push_back(static_cast<unsigned char>(ratio * 255.9)); // b
  }
  bmp::generateBitmapImage(data.data(),
                           freqs_mags.begin()->size(),
                           freqs_mags.size(),
                           imageFile.c_str());
}


// for full windows with even number of points
template<typename T>
std::vector<T> half_hann_window(int const sz) {
  Assert(sz > 0);
  std::vector<T> res;
  res.reserve(sz);
  T const factor = static_cast<T>(M_PI_2) / (2 * sz);
  for (int i=0; i<sz; ++i) {
    res.push_back(std::cos(factor * (1 + 2*i)));
  }
  return res;
}

} // NS
