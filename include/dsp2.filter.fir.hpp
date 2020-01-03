

namespace imajuscule {

struct DescFIRFilter {
    static constexpr int nCoefficientsFadeIn = 0;
    static constexpr bool has_subsampling = false;
    static constexpr bool step_can_error = false;
};

template<typename T, typename Tag>
struct AlgoFIRFilter;

template<typename T, typename Tag>
struct StateFIRFilter {
    
    using Desc = DescFIRFilter;
    using FPT = T;
    using Algo = AlgoFIRFilter<T, Tag>;

    MinSizeRequirement setCoefficients(Algo const & algo, a64::vector<T> v) {
        std::reverse(v.begin(), v.end());
        reset();
        reversed_coefficients.reserve(v.size());
        for(auto const & vv : v) {
            reversed_coefficients.push_back(typename decltype(reversed_coefficients)::value_type(vv));
        }
    
        return {
            static_cast<int>(reversed_coefficients.size()),
            1,
            {}
        };
    }
    
    auto size()  const { return reversed_coefficients.size(); }
    bool isZero() const { return reversed_coefficients.empty(); }
    
    void reset() {
        reversed_coefficients.clear();
    }
    void flushToSilence() {
    }
    
    double getEpsilon(Algo const & algo) const {
        return 2 * std::numeric_limits<T>::epsilon() * reversed_coefficients.size();
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

    struct SetupParam : public Cost {
        void logSubReport(std::ostream & os) const override {
            os << "Brute" << std::endl;
        }
    };
    
    void setup(SetupParam const &) const {}

    constexpr int getLatency() const { return 0; }
    constexpr int getStepPeriod() const { return 1; }

    bool isValid() const {
        return true;
    }
        
    void step(State & s,
              XAndFFTS<T, Tag> const & x_and_ffts,
              Y<T, Tag> & y)
    {
        int constexpr nFrames = 1;
        if(unlikely(s.isZero())) {
            fft::RealSignal_<Tag, T>::zero_n_raw(&y.y[y.progress], nFrames);
            return;
        }
        auto * coeff = s.getReversedCoeffs().data();
        int const sz = static_cast<int>(s.getReversedCoeffs().size());
        for(int i=nFrames-1; i>=0; --i) {
            auto [s1, s2] = x_and_ffts.getSegments(i, sz);
            using E = typename fft::RealSignal_<Tag, FPT>::type::value_type;
            E res;
            dotpr(&x_and_ffts.x[s1.first], coeff, &res, s1.second);
            if(unlikely(s2.second)) {
                E res2;
                dotpr(&x_and_ffts.x[s2.first], coeff+s1.second, &res2, s2.second);
                res += res2;
            }
            y.y[y.progress+nFrames-i-1] += res;
        }
    }
};


}
