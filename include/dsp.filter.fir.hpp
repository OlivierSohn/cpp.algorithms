

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
    
    template<Overlap Mode, typename FFTAlgo>
    MinSizeRequirement getMinSizeRequirement(int const maxVectorSz) const {
        return {
            static_cast<int>(n_coeffs + maxVectorSz - 1), // x block size
            maxVectorSz, // y block size
            {},
            0 // work size
        };
    }
    
    constexpr Latency getImpliedLatency() const {
        // commented out because not constexpr
        //Assert(handlesCoefficients());
        return Latency(0);
    }
    
    int getBiggestScale() const {
        return 1;
    }
    
    template<typename F>
    void forEachUsingSameContext(F f) const {
        f(*this);
    }
    
    int n_coeffs;
};

template<typename T, typename Tag>
struct FIRFilterSimulation {
    
    void setup(FIRSetupParam const & p) {
        stepCost =
        fft::RealSignalCosts<Tag, T>::cost_dotpr(p.n_coeffs) +
        costWriteNConsecutive<T>(1);
    }
    double simuStep(XFFTsCostsFactors const & xFFTCostFactors) {
        return stepCost;
    }
    
    int getBiggestScale() const {
        return 1;
    }
private:
    double stepCost = 0.;
};

template<typename T, typename Tag>
struct Simulation_<FIRSetupParam, T, Tag> {
    using type = FIRFilterSimulation<T, Tag>;
};

template<typename T, typename FFTTag>
struct PartitionAlgo< FIRSetupParam, T, FFTTag> {
    using SetupParam = FIRSetupParam;
    
    static std::optional<SetupParam> run(int n_sources,
                                         int n_channels,
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
