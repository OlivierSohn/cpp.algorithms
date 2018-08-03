

namespace imajuscule
{
  // Always prefer using 'OptimizedFIRFilter' which scales better with the number of coefficients.
  template<typename T>
  struct FIRFilter {
    using FPT = T;
    struct SetupParam {};
    
    void applySetup(SetupParam const &) const {}
    
    void setCoefficients(a64::vector<T> v) {
      past.resize(v.size());
      coefficients = std::move(v);
    }
    
    bool isValid() const {
      return !coefficients.empty();
    }
    
    constexpr int getLatency() const { return 0; }
    
    auto size()  const { return coefficients.size(); }
    bool empty() const { return coefficients.empty(); }

    void clear() {
      coefficients.clear();
      past.resize(0);
    }
    
    T step(T val) {
      if(unlikely(empty())) {
        return {};
      }
      past.feed(val);

      auto res = T{};
      // when coefficients are symmetrical it doesn't matter
      // if we are traversing forward or backward
      // but here we make no assumption:
      auto it_coeff = coefficients.begin();
      // TODO vectorize this (maybe reverse the cycle, so that we can "go forward" in both the input and the coefficients)
      past.for_each_bkwd([&res, &it_coeff](auto val) {
        res += val * *it_coeff;
        ++it_coeff;
      });
      return res;
    }
    
    double getEpsilon() const {
      return 2 * std::numeric_limits<FPT>::epsilon() * coefficients.size();
    }
    
  private:
    a64::vector<T> coefficients;
    cyclic<T> past;
  };
  
  template<typename T>
  static void plotMagnitude(fft::FFTVec<T> const & v) {
    std::vector<T> mags;
    mags.reserve(v.size());
    std::transform(v.begin(), v.end(),
                   std::back_inserter(mags),
                   [](auto v){return abs(v);});
    StringPlot plot(30,1024);
    plot.draw(mags);
    plot.log();
  }
  
  template<typename Container, typename T, typename F>
  auto sample_frequencies(Container & res, bool const even_taps,
                          T const nyquist_freq, F getFreq) {
    
    auto N = res.size();
    auto nyquist = N/2;
    
    T RadPerSample = -M_PI;
    if(even_taps) {
      RadPerSample *= (N - 1.0)/N;
    }
    for(int i=0; i<=nyquist; ++i) {
      auto f = nyquist_freq * i / nyquist;
      T magnitude = getFreq(f);
      auto cplx = magnitude * polar(RadPerSample*i);
      res[i] = cplx;
      if(i && i != nyquist) {
        auto conji = N-i;
        res[conji] = cplx;
      }
    }
    
  }
  
  template<typename T, typename F>
  auto fir_coefficients_by_f_sampling(T nyquist_freq, F getFreq, unsigned int fft_length, unsigned int NumTaps) {
    ScopedLog l("Compute", "FIR coeffs by freq. sampling");
    // according to http://iowahills.com/FIRFiltersByFreqSampling.html
    // with the same number of taps as of fft size (we could try less)
    
    using namespace imajuscule::fft;
    assert(is_power_of_two(fft_length));
    
    a64::vector<complex<T>> res, input;
    res.resize(fft_length);
    input.resize(fft_length);
    
    sample_frequencies(input,
                       (0 == NumTaps%2),
                       nyquist_freq,
                       getFreq);
    
    forward_fft(fft_length, input, res);
    
    auto inv_N = 1. / fft_length;
    
    a64::vector<T> v;
    v.reserve(NumTaps);
    
    // This is where the FIR taps are located in the FFT’s return.
    auto StartT = fft_length/2 - NumTaps/2;
    auto fft_cut_begin = res.begin() + StartT;
    auto fft_cut_end = fft_cut_begin + NumTaps;
    std::transform(fft_cut_begin, fft_cut_end,
                   std::back_inserter(v),
                   [inv_N](auto value) { return value.real()*inv_N; } );
    
    apply_hann_window(v.begin(), v.end());
    
    //plotMagnitude(res);
    return v;
  }
  
  
}
