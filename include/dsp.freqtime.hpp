namespace imajuscule::audio {

/* Throughout this file, you can stride the input signal (parameters named 'windowed_signal_stride').
 If you use a stride > 1 for the input signal, frequency aliasing may occur if your signal has high frequencies in it.
 To avoid this you can either:
 - low-pass your signal before calling these functions,
 - or downsample your signal using 'resampleSinc' and use a value of 1 for 'windowed_signal_stride'
 */

template<typename T>
struct FrequenciesSqMag {
  // indexed by frequency bin:
  std::vector<T> frequencies_sqmag;
  
  int get_fft_length() const {
    return fft_length;
  }
  void set_fft_length(int l) {
    fft_length = l;
  }
  double bin_index_to_Hz(int const samplingRate) const {
    Assert(fft_length);
    return samplingRate / static_cast<double>(fft_length);
  }

private:
  int fft_length; // to convert frequency from "bin index" representation to Herz representation
};

template<typename ITER>
void apply_window(ITER it,
                  ITER const end,
                  int const signal_stride,
                  std::vector<typename ITER::value_type> const & half_window,
                  a64::vector<typename ITER::value_type> & v) {
  v.clear();
  Assert(2*half_window.size() <= v.capacity());
  auto wEnd = half_window.end();
  auto wBegin = half_window.begin();
  auto wIt = wEnd - 1;
  bool first_half = true;
  for(; it < end; it += signal_stride) {
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
  Assert(wIt == wEnd);
  Assert(it >= end);
}


template<typename ITER>
void apply_rectangular_window(ITER it,
                              ITER const end,
                              int const signal_stride,
                              a64::vector<typename ITER::value_type> & v) {
  v.clear();
  for(; it < end; it += signal_stride) {
    v.push_back(*it);
  }
}

inline int get_fft_length_for(int const full_window_sz,
                              int const zero_padding_factor) {
  return ceil_power_of_two(full_window_sz * zero_padding_factor);
}

template<typename Tag, typename ITER>
void findFrequenciesSqMag(ITER it,
                          ITER const end,
                          int const windowed_signal_stride,
                          std::vector<typename ITER::value_type> const & half_window,
                          int const zero_padding_factor,
                          typename fft::RealSignal_<Tag, typename ITER::value_type>::type & signal,
                          typename fft::RealFBins_<Tag, typename ITER::value_type, a64::Alloc>::type & result,
                          FrequenciesSqMag<typename ITER::value_type> & frequencies_sqmag,
                          std::function<void(a64::vector<typename ITER::value_type>)> fDebugWindowedSignal = {}) {
  {
    using VAL = typename ITER::value_type;
    using Algo = fft::Algo_<Tag, VAL>;
    
    int const numSamplesUnstrided = static_cast<int>(std::distance(it, end));
    int const numSamplesStrided = 1 + (numSamplesUnstrided-1) / windowed_signal_stride;
    Assert(numSamplesStrided == static_cast<int>(2 * half_window.size()));
    int const fft_length = get_fft_length_for(numSamplesStrided, zero_padding_factor);

    using Contexts = fft::Contexts_<Tag, VAL>;
    Algo fft(Contexts::getInstance().getBySize(fft_length));
    
    signal.clear();
    reserve_no_shrink(signal,
                      fft_length);
    
    // apply the window...
    apply_window(it, end, windowed_signal_stride, half_window, signal);
    
    // ... and pad
    signal.resize(fft_length, VAL{0});
    
    if (fDebugWindowedSignal) {
      fDebugWindowedSignal(signal);
    }
    
    result.resize(fft_length);
    fft.forward(signal.begin(), result.data(), fft_length);
    fft::unwrap_frequencies_sqmag<Tag>(result, fft_length, frequencies_sqmag.frequencies_sqmag);
    frequencies_sqmag.set_fft_length(fft_length);
    
    const VAL factor = 1. / (Algo::scale * Algo::scale * fft_length * fft_length);
    for (auto & v : frequencies_sqmag.frequencies_sqmag) {
      v *= factor;
    }
  }
}

/** Returns the frequencies of a windowed signal.
 *
 * For better efficiency, the signal length should be a power of 2
 * @param input_stride stride for input
 * @param half_window is the right window half, starting at 1.0 (hence the full window has 2 times the 1.0 value)
 * @param frequencies_sqmag is the sqared magnitude of the frequencies from 0 to nyquist freq included.
 */
template<typename ITER>
void findFrequenciesSqMagSlow(ITER it,
                              ITER const end,
                              int const windowed_signal_stride,
                              std::vector<typename ITER::value_type> const & half_window,
                              int const zero_padding_factor,
                              FrequenciesSqMag<typename ITER::value_type> & frequencies_sqmag,
                              std::function<void(a64::vector<typename ITER::value_type>)> fDebugWindowedSignal = {}) {
  using VAL = typename ITER::value_type;
  using Tag = fft::Fastest;
  
  using RealFBins = fft::RealFBins_<Tag, VAL, a64::Alloc>;
  
  typename RealFBins::type result;
  typename fft::RealSignal_<Tag, VAL>::type signal;
  
  findFrequenciesSqMag<Tag>(it,
                            end,
                            windowed_signal_stride,
                            half_window,
                            zero_padding_factor,
                            signal,
                            result,
                            frequencies_sqmag,
                            fDebugWindowedSignal);
}

template<typename T>
struct FreqMag {
  T freq; // unit: Herz
  T mag_db;
  
  bool operator == (FreqMag const & o) const {
    return std::make_tuple(freq, mag_db) ==
    std::make_tuple(o.freq, o.mag_db);
  }
};

template<typename T>
std::ostream & operator << (std::ostream& os, FreqMag<T> const & f) {
  os << "FreqMag(f= " << f.freq << " db= " << f.mag_db << ")";
  return os;
}

template<typename T>
struct QuadraticInterpolation {
  QuadraticInterpolation(T const alpha, T const beta, T const gamma) {
    // extractMaxFreqSqMag enforces:
    // betaAmplitude >= alphaAmplitude && betaAmplitude > gammaAmplitude
    // but then the amplitude -> dB transformation can turn strict inequalities into weak ones.
    Assert(beta >= alpha && beta >= gamma);
    
    T denom = alpha - 2.*beta + gamma;
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


// 'transform_amplitude' must be such that quadratic interpolation "works well" for peak finding
template<typename T, typename Allocator, typename TransformAmplitude, typename F>
void foreachLocalMaxFreqsMags(std::vector<T, Allocator> const & freqs_sqmag,
                              TransformAmplitude transform_amplitude,
                              F f) {
  Assert(!freqs_sqmag.empty());
  
  for (int i=0, sz=static_cast<int>(freqs_sqmag.size()); i<sz; ++i) {
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
    std::optional<T> const beta = transform_amplitude(freqs_sqmag[i]);
    if (!beta) {
      continue;
    }
    if (unlikely(i == 0 || i == sz-1)) {
      // we cannot interpolate, we are at a boundary
      f(i,
        *beta);
    } else {
      // use quadratic interpolation
      auto alpha = transform_amplitude(freqs_sqmag[i-1]);
      auto gamma = transform_amplitude(freqs_sqmag[i+1]);
      if (!alpha || !gamma) {
        // we cannot interpolate
        f(i,
          *beta);
      }
      auto itp = QuadraticInterpolation(*alpha,
                                        *beta,
                                        *gamma);
      f(i + itp.getP(),
        itp.getMag());
    }
  }
}

template<typename T, typename TransformAmplitude>
void extractLocalMaxFreqsMags(double const samplingRate,
                              FrequenciesSqMag<T> const & freqs_sqmag,
                              TransformAmplitude transform_amplitude,
                              std::vector<FreqMag<T>> & res) {
  res.clear();
  reserve_no_shrink(res,
                    freqs_sqmag.frequencies_sqmag.size() / 2); // pessimistically
  
  double const bin_index_to_Hz = freqs_sqmag.bin_index_to_Hz(samplingRate);
  foreachLocalMaxFreqsMags(freqs_sqmag.frequencies_sqmag,
                           transform_amplitude,
                           [&res, bin_index_to_Hz](T const frequency, T const magnitude) {
    res.push_back({
      bin_index_to_Hz * frequency,
      magnitude
    });
  });
}

/* Converts a square magnitude (homogenous to power)
 * to decibels
 */
template<typename T>
struct SqMagToDb {
  std::optional<T> operator() (T const val){
    Assert(val >= 0);
    if (val == 0.) {
      return {};
    }
    return 10. * std::log10(val);
  }
};
template<typename T>
struct DbToSqMag {
  T operator() (T const val){
    return std::pow(10., val * 0.1);
  }
};

template<typename T>
struct MagToDb {
  std::optional<T> operator() (T const val){
    Assert(val >= 0);
    if (val == 0.) {
      return {};
    }
    return 20. * std::log10(val);
  }
};
template<typename T>
struct DbToMag {
  T operator() (T const val){
    // using
    // sqrt(x) = x^(0.5)
    // and
    // (a^b)^c = a^(b*c)
    return std::pow(10., val * 0.05);
  }
};


template<typename Iter>
void findLocalMaxMagFrequenciesSlow(Iter it,
                                    Iter const end,
                                    double const signal_sampling_rate,
                                    int const windowed_signal_stride,
                                    std::vector<typename Iter::value_type> const & half_window,
                                    int window_center_stride,
                                    int const zero_padding_factor,
                                    std::vector<std::vector<FreqMag<typename Iter::value_type>>> & freqs_sqmags) {
  using T = typename Iter::value_type;
  
  int const szWindow = 2 * half_window.size();
  int const windowSampleSpan = (szWindow - 1) * windowed_signal_stride + 1;
  int const numSamples = std::distance(it, end);
  int const num_ffts = 1 + (numSamples - windowSampleSpan) / window_center_stride;
  
  auto limit = it + windowSampleSpan;
  
  freqs_sqmags.clear();
  freqs_sqmags.resize(num_ffts);
  auto itRes = freqs_sqmags.begin();
  
  FrequenciesSqMag<T> freqs_sqmag;
  int control = 0;
  for (; limit <= end; it += window_center_stride, limit += window_center_stride, ++itRes) {
    ++control;
    
    findFrequenciesSqMagSlow(it, limit, windowed_signal_stride, half_window, zero_padding_factor, freqs_sqmag
#if 0
                             , [control, it, limit](a64::vector<T> const & windowed_signal) {
      std::string ctrl;
      if (control < 100) {
        ctrl += "0";
      }
      if (control < 10) {
        ctrl += "0";
      }
      ctrl += std::to_string(control);
      
      std::vector<T> signal(it, limit);
      
      bool zero = true;
      for (auto v : signal) {
        if (v) {
          zero = false;
        }
      }
      for (auto v : windowed_signal) {
        if (v) {
          zero = false;
        }
      }
      if (zero) {
        return;
      }
      write_wav("/Users/Olivier/",
                ctrl + ".signal.findLocalMaxMagFrequencies.wav",
                std::vector<std::vector<T>>{signal},
                100);
      write_wav("/Users/Olivier/",
                ctrl + ".windowedsignal.findLocalMaxMagFrequencies.wav",
                std::vector<a64::vector<T>>{windowed_signal},
                100);
    }
#endif
                             );
    extractLocalMaxFreqsMags(signal_sampling_rate / windowed_signal_stride,
                             freqs_sqmag,
                             SqMagToDb<T>(),
                             *itRes);
  }
  Assert(control == num_ffts);
}

template<typename Iter>
void findAllMagFrequenciesSlow(Iter it,
                               Iter const end,
                               int const windowed_signal_stride,
                               std::vector<typename Iter::value_type> const & half_window,
                               int window_center_stride,
                               int const zero_padding_factor,
                               std::vector<FrequenciesSqMag<typename Iter::value_type>> & freqs_mags) {
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
    findFrequenciesSqMagSlow(it, limit, windowed_signal_stride, half_window, zero_padding_factor, *itRes);
  }
  Assert(freqs_mags.size() == num_ffts);
}

template<typename Iter>
void drawSpectrumSlow(Iter it,
                      Iter const end,
                      int const windowed_signal_stride,
                      std::vector<typename Iter::value_type> const & half_window,
                      int window_center_stride,
                      int const zero_padding_factor,
                      std::string const& imageFile) {
  using T = typename Iter::value_type;
  std::vector<FrequenciesSqMag<T>> freqs_mags;
  
  findAllMagFrequenciesSlow(it,
                            end,
                            windowed_signal_stride,
                            half_window,
                            window_center_stride,
                            zero_padding_factor,
                            freqs_mags);
  
  std::vector<T> all_freqs_mags;
  reserve_no_shrink(all_freqs_mags,
                    freqs_mags.size() * freqs_mags.begin()->frequencies_sqmag.size());
  for (int i=0, sz = freqs_mags.begin()->frequencies_sqmag.size(); i<sz; ++i) {
    for (auto const & v : freqs_mags) {
      all_freqs_mags.push_back(v.frequencies_sqmag[i]);
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
                           freqs_mags.begin()->frequencies_sqmag.size(),
                           freqs_mags.size(),
                           imageFile.c_str());
}

template<typename T>
struct DeducedNote {
  static_assert(std::is_floating_point_v<T>);
  
  T frequency;
  T amplitude;
  
  // frame units : window stride
  int startFrame;
  int endFrame;
};

template<typename T>
std::ostream & operator << (std::ostream & os, DeducedNote<T> const & n) {
  os << "DeducedNote( f " << std::setw(10) << n.frequency
  << ", amp " << std::setw(10) << n.amplitude
  << ", range " << std::setw(5) << n.startFrame
  << " " << std::setw(5) << n.endFrame << ")";
  return os;
}


template<typename T>
void drawDeducedNotes(std::vector<DeducedNote<T>> const & notes,
                      double const lowest_detectable_frequency,
                      std::string const& imageFile) {
  Optional<T> max_mag, min_mag;
  for (auto const & val : notes) {
    if (!max_mag || val.amplitude > *max_mag) {
      max_mag = val.amplitude;
    }
    if (!min_mag || val.amplitude < *min_mag) {
      min_mag = val.amplitude;
    }
  }
  
  Assert(max_mag);
  Assert(min_mag);
  
  auto frequencyToImageHeight = [](T freq) {
    return static_cast<int>(freq);
  };
  
  Optional<T> max_f = lowest_detectable_frequency;
  Optional<int> max_frame, min_frame;
  
  std::multimap<T,DeducedNote<T>> byAmplitude, byFreq;
  std::multimap<int,DeducedNote<T>> byDuration, byStart;
  for (auto const & val : notes) {
    byFreq.emplace(val.frequency, val);
    byAmplitude.emplace(val.amplitude, val);
    byDuration.emplace(1 + val.endFrame - val.startFrame, val);
    byStart.emplace(val.startFrame, val);
    if (!max_f || val.frequency > *max_f) {
      max_f = val.frequency;
    }
    if (!max_frame || val.endFrame > *max_frame) {
      max_frame = val.endFrame;
    }
    if (!min_frame || val.startFrame < *min_frame) {
      min_frame = val.startFrame;
    }
  }
  Assert(max_f);
  Assert(max_frame);
  Assert(min_frame);
  
  T freq_span = *max_f;
  T mag_span = *max_mag - *min_mag;
  
  std::vector<std::vector<T>> amplitudes;
  
  amplitudes.resize(1 + *max_frame - *min_frame);
  for (auto & v: amplitudes) {
    v.resize(static_cast<int>(freq_span) + 1,
             0.);
  }
  
  for (auto const & n : notes) {
    for (int i = n.startFrame; i <= n.endFrame; ++i) {
      auto & a = amplitudes[i-*min_frame][frequencyToImageHeight(n.frequency)];
      a = std::max(a,
                   mag_span? (1.0 - (*max_mag - n.amplitude) / mag_span) : 0.5);
    }
  }
  
  std::vector<T> all_amplitudes;
  std::vector<bool> red;
  all_amplitudes.reserve(amplitudes.size() * amplitudes.begin()->size());
  red.reserve(amplitudes.size() * amplitudes.begin()->size());
  for (int i=0, sz = amplitudes.begin()->size(); i<sz; ++i) {
    for (auto const & v : amplitudes) {
      if (i==frequencyToImageHeight(lowest_detectable_frequency)) {
        red.push_back(true);
      } else {
        red.push_back(false);
      }
      all_amplitudes.push_back(v[i]);
    }
  }
  
  
  std::vector<unsigned char> data;
  data.reserve(3 * all_amplitudes.size());
  for (auto const & ratio : all_amplitudes) {
    Assert(ratio >= 0.f);
    Assert(ratio <= 1.f);
    data.push_back(static_cast<unsigned char>(ratio * 255.9)); // r
    data.push_back(static_cast<unsigned char>(ratio * 255.9)); // g
    data.push_back(static_cast<unsigned char>(ratio * 255.9)); // b
  }
  for (int i=0, sz=static_cast<int>(red.size()); i<sz; ++i) {
    if (red[i]) {
      data[3*i+2] = 255;
    }
  }
  
  bmp::generateBitmapImage(data.data(),
                           amplitudes.begin()->size(),
                           amplitudes.size(),
                           imageFile.c_str());
  
  std::cout << "deduced notes by amplitude:" << std::endl;
  for (auto const & v : byAmplitude) {
    std::cout << v.second << std::endl;
  }  std::cout << "deduced notes by freq:" << std::endl;
  for (auto const & v : byFreq) {
    std::cout << v.second << std::endl;
  }
  std::cout << "deduced notes by duration:" << std::endl;
  for (auto const & v : byDuration) {
    std::cout << v.second << std::endl;
  }
#if 0 // takes a long time
  std::cout << "deduced notes by start time:" << std::endl;
  for (auto const & v : byStart) {
    std::cout << "|";
    int count = 0;
    {
      int nBlanks = v.second.startFrame - *min_frame;
      int n10Blanks = nBlanks / 10;
      nBlanks -= 10 * n10Blanks;
      std::cout << std::setw(5) << 10 * n10Blanks;
      std::cout << "|";
      for (int i = 0; i < nBlanks; ++i) {
        std::cout << " ";
        ++count;
      }
    }
    for (int i=v.second.startFrame; i <= v.second.endFrame; ++i) {
      std::cout << "-";
      ++count;
    }
    
    {
      int nBlanks = *max_frame - v.second.endFrame;
      int n10Blanks = nBlanks / 10;
      nBlanks -= 10 * n10Blanks;
      for (int i = 0; i < nBlanks; ++i) {
        std::cout << " ";
        ++count;
      }
      std::cout << "|";
      std::cout << std::setw(5) << 10 * n10Blanks;
    }
    std::cout << "|";
    int pad = 50-count;
    if (pad > 0) {
      std::cout << std::string(pad, ' ');
    }
    std::cout << v.second << std::endl;
  }
#endif
}


/**
 @param frequencyEpsilon : if abs(ln(f1)-ln(f2)) < frequencyEpsilon, f1 and f2 are considered to be equal.
 Note that the natural log of 2 frequencies that are one half tone apart is ln(2) / 12 = 0.05776226504
 */
template<typename Iter>
auto deduceNotesSlow(Iter it,
                     Iter const end,
                     double const signal_sample_rate,
                     int const windowed_signal_stride,
                     std::vector<typename Iter::value_type> const & half_window,
                     int window_center_stride,
                     int const zero_padding_factor,
                     typename Iter::value_type const frequencyEpsilon)
-> std::vector<DeducedNote<typename Iter::value_type>> {
  using T = typename Iter::value_type;
  std::vector<std::vector<FreqMag<T>>> freqs_mags;
  
  findLocalMaxMagFrequenciesSlow(it,
                                 end,
                                 signal_sample_rate,
                                 windowed_signal_stride,
                                 half_window,
                                 window_center_stride,
                                 zero_padding_factor,
                                 freqs_mags);
  
  std::vector<DeducedNote<T>> res;
  res.reserve(100);
  
  struct AlmostFrequency {
    AlmostFrequency(T t, T epsilon)
    : val(t)
    , epsilon(epsilon) {
    }
    
    bool operator == (AlmostFrequency const & o) const {
      Assert(epsilon == o.epsilon);
      return std::abs(val - o.val) <= epsilon;
    }
    bool operator < (AlmostFrequency const & o) const {
      Assert(epsilon == o.epsilon);
      return val < o.val && !operator == (o);
    }
    bool operator > (AlmostFrequency const & o) const {
      Assert(epsilon == o.epsilon);
      return val > o.val && !operator == (o);
    }
    
  private:
    T val;
    T epsilon;
  };
  
  auto almostFreq = [frequencyEpsilon](T f) { return AlmostFrequency{f, frequencyEpsilon}; };
  
  struct History {
    History(int frame, FreqMag<T> fm)
    : startFrame(frame)
    , endFrame(frame)
    {
      fms.reserve(20);
      fms.push_back(fm);
    }
    
    DeducedNote<T> deduceNote() const {
      DeducedNote<T> r;
      r.startFrame = startFrame;
      r.endFrame = endFrame;
      
      r.frequency = 0.;
      r.amplitude = 0.;
      for (auto const & v : fms) {
        r.frequency += v.freq;
        r.amplitude += DbToMag<double>()(v.mag_db);
      }
      r.frequency /= fms.size();
      r.amplitude /= fms.size();
      return r;
    }
    
    int startFrame, endFrame;
    std::vector<FreqMag<T>> fms;
  };
  
  std::map<AlmostFrequency, History> cur;
  int frame = 0;
  for (int sz=freqs_mags.size(); frame < sz; ++frame) {
    for (FreqMag<T> const & fm : freqs_mags[frame]) {
      auto const f = almostFreq(std::log(fm.freq));
      auto it = cur.find(f);
      if (it == cur.end()) {
        cur.emplace(f,
                    History(frame, fm));
      } else if (it->second.endFrame == frame) {
        // can happen when frequency epsilon is significant
        it->second.fms.push_back(fm);
      } else {
        Assert(it->second.endFrame == frame-1);
        it->second.endFrame = frame;
        it->second.fms.push_back(fm);
      }
    }
    // end some notes if needed
    for (auto it = cur.begin(); it != cur.end();) {
      Assert(it->second.endFrame <= frame);
      if (it->second.endFrame < frame) {
        Assert(it->second.endFrame == frame - 1);
        res.push_back(it->second.deduceNote());
        it = cur.erase(it);
      } else {
        ++it;
      }
    }
  }
  
  for (auto it = cur.begin(); it != cur.end(); ++it) {
    Assert(it->second.endFrame == frame-1);
    res.push_back(it->second.deduceNote());
  }
  return res;
}


// for full windows with even number of points
template<typename T>
std::vector<T> half_hann_window(int const half_sz) {
  Assert(half_sz > 0);
  std::vector<T> res;
  res.reserve(half_sz);
  T const factor = static_cast<T>(M_PI_2) / (2 * half_sz);
  for (int i=0; i<half_sz; ++i) {
    res.push_back(std::cos(factor * (1 + 2*i)));
  }
  return res;
}
template<typename T>
void half_gaussian_window(int sigma_factor,
                          int const half_sz,
                          std::vector<T> & res) {
  Assert(half_sz > 0);
  res.clear();
  reserve_no_shrink(res,
                    half_sz);
  T const maxT = sigma_factor;
  T const increment = maxT / half_sz;
  
  // with sigma = 1
  
  for (int i=0; i<half_sz; ++i) {
    T const t = increment * (i + 0.5);
    res.push_back(std::exp(-(t*t) * 0.5));
  }
}
template<typename T>
void half_rectangular_window(int const half_sz,
                             std::vector<T> & res) {
  Assert(half_sz > 0);
  res.clear();
  reserve_no_shrink(res,
                    half_sz);

  for (int i=0; i<half_sz; ++i) {
    res.push_back(1.);
  }
}

template<typename T>
void normalize_window(std::vector<T> & w) {
  if (w.empty()) {
    return;
  }
  T avg{};
  for (auto const & v : w) {
    // According to tests (see TEST(FFT, unwrap_freq_sqmag)),
    // to find the right _magnitude_ in the spectrum, we can use the average on the first order:
    avg += v;
    
    // see also https://community.sw.siemens.com/s/article/window-correction-factors
  }
  avg /= w.size();
  Assert(avg);
  T const inv_avg = 1. / avg;
  std::transform(w.begin(),
                 w.end(),
                 w.begin(),
                 [inv_avg](T const & val) { return val * inv_avg; });
}



/*
 when 2 signals are correlated, an equal-gain crossfade must be used.
 when 2 signals are not correlated, an equal-power crossfade must be used
 */

enum class EqualGainCrossFade {
  Linear,    // first derivative is discontinuous
  Sinusoidal // all derivatives are continuous
};

enum class EqualPowerCrossFade {
  
};

template<typename T>
struct XFadeValues {
  T new_signal_mult;
  T old_signal_mult;
};

// The equal-gain crossfade has the property : f(x) = 1 - f(1-x)
// so we store the fist half, and deduce the second half using this property.
template<typename T>
struct EqualGainXFade {

  void reserve(std::size_t sz) {
    reserve_no_shrink(first_half,
                      sz/2);
  }

  // @param fullsize : the number of frames where the 2 signals are mixed with non-zero gains
  //                   (must be even)
  void set(int const full_size,
           EqualGainCrossFade const t) {
    // we could support odd-length crossfades but the code would be more complicated
    Assert(0 == full_size % 2);
    if (type && *type == t && full_size == 2*half_size) {
      return;
    }
    type = t;
    
    half_size = full_size / 2;

    first_half.clear();
    reserve_no_shrink(first_half,
                      half_size);

    for (int i=1; i <= half_size; ++i) {
      // when full_size = 2, with a linear fade we want:
      // 1/3 2/3
      T const x = i / static_cast<T>(full_size+1);
      T y;
      switch(t) {
        case EqualGainCrossFade::Linear:
          y = x;
          break;
        case EqualGainCrossFade::Sinusoidal:
          y = 0.5 * (1. + std::sin(M_PI * x - M_PI_2));
          break;
      }
      first_half.push_back(y);
    }
  }

  XFadeValues<T> get(int const progress) const {
    if (progress <= 0) {
      return {0., 1.};
    } else if (progress <= half_size) {
      auto const v = first_half[progress-1];
      return {v, 1. - v};
    } else if (progress <= 2 * half_size) {
      auto const v = first_half[(2*half_size)-progress];
      return {1. - v, v};
    } else {
      return {1., 0.};
    }
  }
private:
  int half_size;
  std::vector<T> first_half;
  std::optional<EqualGainCrossFade> type;
};

} // NS
