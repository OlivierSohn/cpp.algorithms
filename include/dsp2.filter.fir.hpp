

namespace imajuscule {

struct DescFIRFilter {
    static constexpr int nCoefficientsFadeIn = 0;
    static constexpr bool has_subsampling = false;
    static constexpr bool step_can_error = false;
};

template<typename T, typename FFTTag>
struct AlgoFIRFilter;

template<typename T, typename FFTTag>
struct StateFIRFilter {
    using Algo = AlgoFIRFilter<T, FFTTag>;
    using Desc = DescFIRFilter;
    
    using FPT = T;
    using Tag = FFTTag;

    MinSizeRequirement setCoefficients(Algo const & algo,
                                       a64::vector<T> v)
    {
        reset();

        std::reverse(v.begin(), v.end());
        reversed_coefficients.reserve(v.size());
        for(auto const & vv : v) {
            reversed_coefficients.push_back(typename decltype(reversed_coefficients)::value_type(vv));
        }
    
        return {
            static_cast<int>(reversed_coefficients.size()), // x block size
            1, // y block size
            0, // anticipated y writes
            {}
        };
    }
        
    template<typename F>
    void onContextFronteer(F f) {
    }

    auto size()  const { return reversed_coefficients.size(); }
    bool isZero() const { return reversed_coefficients.empty(); }
    
    void reset() {
        // if needed, to avoid singularities, we could put a single zero coefficient here
        reversed_coefficients.clear();
    }
    
    double getEpsilon(Algo const & algo) const {
        return 2 * std::numeric_limits<T>::epsilon() * reversed_coefficients.size();
    }
    void logComputeState(Algo const & algo, std::ostream & os) const {
        os << "Brute [_/" << reversed_coefficients.size() << "]" << std::endl;
    }
    using RealSignal = typename fft::RealSignal_<Tag, T>::type;
    RealSignal const & getReversedCoeffs() const { return reversed_coefficients; }
private:
    RealSignal reversed_coefficients;
};


template<typename T, typename FFTTag>
struct AlgoFIRFilter {
    using Desc = DescFIRFilter;
    using State = StateFIRFilter<T, FFTTag>;
    using FPT = T;
    using Tag = FFTTag;
    
    static constexpr auto dotpr = fft::RealSignal_<Tag, FPT>::dotpr;
    
    using SetupParam = FIRSetupParam;

    void setup(SetupParam const &) const {}

    constexpr Latency getLatency() const {
        // commented out because not constexpr
        //Assert(handlesCoefficients());
        return Latency(0);
    }
    
    constexpr bool handlesCoefficients() const {
        return true;
    }
    bool isValid() const {
        return true;
    }

    void dephaseSteps(State &s,
                      int n_steps) const {
    }

    void step(State & s,
              XAndFFTS<T, Tag> const & x_and_ffts,
              Y<T, Tag> & y) const
    {
        auto * coeff = s.getReversedCoeffs().data();

        int const sz = static_cast<int>(s.getReversedCoeffs().size());

        auto [s1, s2] = x_and_ffts.getSegments(0, sz);

        using E = typename fft::RealSignal_<Tag, FPT>::type::value_type;
        E res; // no need to zero-initialize, dotpr will write it.

        // s1.second may be 0
        dotpr(&x_and_ffts.x[s1.first], coeff, &res, s1.second);

        if(unlikely(s2.second)) {
            E res2;
            dotpr(&x_and_ffts.x[s2.first], coeff+s1.second, &res2, s2.second);
            res += res2;
        }

        y.y[y.uProgress] += res;
    }

    void flushToSilence(State & s) const {
    }
};

template<typename T, typename FFTTag>
struct corresponding_legacy_dsp<AlgoFIRFilter<T, FFTTag>> {
    // FFTTag is ignored, legacy FIRFilter uses the fastest one available
    using type = FIRFilter<T>;
};

}
