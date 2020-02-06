

namespace imajuscule
{
struct FIRSetupParam : public Cost {
    FIRSetupParam(int n_coeffs)
    : n_coeffs(n_coeffs)
    {}
    
    void logSubReport(std::ostream & os) const override {
        os << "Brute " << n_coeffs << " coeffs" << std::endl;
    }

    constexpr bool isValid() const {
        return true;
    }
    constexpr bool handlesCoefficients() const {
        return n_coeffs > 0;
    }
    
    void adjustWork(int targetNCoeffs) {
        Assert(n_coeffs >= targetNCoeffs);
        n_coeffs = std::min(n_coeffs, targetNCoeffs);

        if(!handlesCoefficients()) {
            setCost(0.);
        }
    }
    
    MinSizeRequirement getMinSizeRequirement() const {
        return {
            static_cast<int>(n_coeffs), // x block size
            1, // y block size
            0, // anticipated y writes
            {},
            0 // work size
        };
    }
    
    constexpr Latency getImpliedLatency() const {
        // commented out because not constexpr
        //Assert(handlesCoefficients());
        return Latency(0);
    }
    
    template<typename F>
    void forEachUsingSameContext(F f) const {
        f(*this);
    }
    
    int n_coeffs;
};

  /*
   Brute force filtering (no fft is used).
   */
  template<typename T, template<typename> typename Allocator>
  struct FIRFilter {
    using FPT = T;
    static constexpr int nComputePhaseable = 0;
    static constexpr int nCoefficientsFadeIn = 0;
    static constexpr bool has_subsampling = false;
    static constexpr bool step_can_error = false;

    using SetupParam = FIRSetupParam;

    void logComputeState(std::ostream & os) const {
        os << "Brute [" << past.getIndex() << "/" << reversed_coefficients.size() << "]" << std::endl;
    }
    static constexpr auto dotpr = fft::RealSignal_<fft::Fastest, FPT>::dotpr;
    
    void setup(SetupParam const &) const {}

      std::array<int, nComputePhaseable> getComputePeriodicities() const {
          return {};
      }
      // in [0, getComputePeriodicity())
      std::array<int, nComputePhaseable> getComputeProgresses() const {
          return {};
      }
      void setComputeProgresses(std::array<int, nComputePhaseable> const & progresses) {
      }
      
      void setCoefficients(a64::vector<T> v) {
          past.resize(v.size());
          std::reverse(v.begin(), v.end());
          reversed_coefficients = std::move(v);
      }
    
    bool handlesCoefficients() const {
        return true;
    }
    constexpr bool isValid() const {
      return true;
    }
    
    constexpr Latency getLatency() const {
        Assert(handlesCoefficients());
        return Latency(0);
    }
    
    auto size()  const { return reversed_coefficients.size(); }
    bool isZero() const { return reversed_coefficients.empty(); }

    void reset() {
      reversed_coefficients.clear();
      past.resize(0);
    }
    void flushToSilence() {
      past.reset();
    }
    
    T step(T val) {
      if(unlikely(isZero())) {
        return {};
      }
      return doStep(val);
    }
  private:
      T doStep(T val) {
          past.feed(val);
          
          auto [s1, s2] = past.forward_traversal();
          auto * coeff = reversed_coefficients.data();
          T res1, res2;
          dotpr(s1.first, coeff,           &res1, s1.second);
          dotpr(s2.first, coeff+s1.second, &res2, s2.second);
          return res1 + res2;
      }
  public:
      template<typename FPT2>
      void stepAssignVectorized(FPT2 const * const input_buffer,
                                FPT2 * output_buffer,
                                int nSamples)
      {
          if(unlikely(isZero())) {
              // zero output_buffer
              fft::RealSignal_<fft::Fastest, FPT2>::zero_n_raw(output_buffer, nSamples);
              return;
          }
          for(int i=0; i<nSamples; ++i) {
              output_buffer[i] = doStep(input_buffer[i]);
          }
      }
      template<typename FPT2>
      void stepAddVectorized(FPT2 const * const input_buffer,
                                FPT2 * output_buffer,
                                int nSamples)
      {
          if(unlikely(isZero())) {
            return;
          }
          for(int i=0; i<nSamples; ++i) {
              output_buffer[i] += doStep(input_buffer[i]);
          }
      }
      template<typename FPT2>
      void stepAddInputZeroVectorized(FPT2 * output_buffer,
                                      int nSamples)
      {
          if(unlikely(isZero())) {
            return;
          }
          for(int i=0; i<nSamples; ++i) {
              output_buffer[i] += doStep({});
          }
      }

    double getEpsilon() const {
      return 2 * std::numeric_limits<FPT>::epsilon() * reversed_coefficients.size();
    }
    
  private:
    std::vector<T, Allocator<T>> reversed_coefficients;
    cyclic<T> past;
  };

template<typename T, typename FFTTag>
struct PartitionAlgo< FIRSetupParam, T, FFTTag> {
    using SetupParam = FIRSetupParam;
    
    static std::optional<SetupParam> run(int n_channels,
                                         int n_audio_channels,
                                         int n_audio_frames_per_cb,
                                         int total_response_size,
                                         int n_scales,
                                         double frame_rate,
                                         std::ostream & os) {
        os << "Optimization of FIRFilter" << std::endl;
        IndentingOStreambuf i(os);

        // there is no variable to optimize with FIRFilter:
        SetupParam p{total_response_size};
        p.setCost(0.);
        return p;
    }
};


  template<typename T>
  static void plotMagnitude(fft::FFTVec<T> const & v) {
    std::vector<T> mags;
    mags.resize(v.size());
    std::transform(v.begin(), v.end(),
                   mags.begin(),
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
        
    // This is where the FIR taps are located in the FFT’s return.
    auto StartT = fft_length/2 - NumTaps/2;
    auto fft_cut_begin = res.begin() + StartT;
    auto fft_cut_end = fft_cut_begin + NumTaps;
    a64::vector<T> v;
    v.resize(NumTaps);
    std::transform(fft_cut_begin, fft_cut_end,
                   v.begin(),
                   [inv_N](auto value) { return value.real()*inv_N; } );
    
    apply_hann_window(v.begin(), v.end());
    
    //plotMagnitude(res);
    return v;
  }
  
  
}
