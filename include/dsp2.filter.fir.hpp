namespace imajuscule::audio {

struct DescFIRFilter {
    static constexpr int nCoefficientsFadeIn = 0;
    static constexpr bool has_subsampling = false;
    static constexpr bool step_can_error = false;
};

template<typename T, template<typename> typename Allocator, typename FFTTag>
struct AlgoFIRFilter;

template<typename T, template<typename> typename Allocator, typename FFTTag>
struct StateFIRFilter {
    using Algo = AlgoFIRFilter<T, Allocator, FFTTag>;
    using Desc = DescFIRFilter;
    
    using FPT = T;
    using Tag = FFTTag;

    static int getAllocationSz_SetCoefficients(typename Algo::SetupParam const & p) {
        return 0;
    }

    void setCoefficients(Algo const & algo,
                         a64::vector<T> v)
    {
        reset();

        if(v.size() != algo.n_coeffs) {
            throw std::logic_error("wrong count of coeffs");
        }
        std::reverse(v.begin(), v.end());
        reversed_coefficients.reserve(v.size());
        for(auto const & vv : v) {
            reversed_coefficients.push_back(typename decltype(reversed_coefficients)::value_type(vv));
        }
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
    using RealSignalOutput = typename fft::RealSignal_<Tag, T>::outputType;
    RealSignalOutput const & getReversedCoeffs() const { return reversed_coefficients; }
private:
    RealSignalOutput reversed_coefficients;
};


template<typename T, template<typename> typename Alloc, typename FFTTag>
struct AlgoFIRFilter {
    template<typename TT>
    using Allocator = Alloc<TT>;
    
    using Desc = DescFIRFilter;
    using State = StateFIRFilter<T, Alloc, FFTTag>;
    using FPT = T;
    using Tag = FFTTag;
    
    static constexpr auto dotpr = fft::RealSignal_<Tag, FPT>::dotpr;
    
    using SetupParam = FIRSetupParam;

    void setup(SetupParam const &p) {
        n_coeffs = p.n_coeffs;
    }

    constexpr Latency getLatency() const {
        // commented out because not constexpr
        //Assert(handlesCoefficients());
        return Latency(0);
    }

    int getBiggestScale() const {
        return 1;
    }
    
    bool handlesCoefficients() const {
        return n_coeffs > 0;
    }
    bool isValid() const {
        return n_coeffs >= 0;
    }

    void dephaseStep(State &s,
                     int x_progress) const {
    }
    
    template<template<typename> typename Allocator2, typename WorkData>
    void stepVectorized(State & state,
                        XAndFFTS<T, Allocator2, Tag> const & x_and_ffts,
                        Y<T, Tag> & y,
                        WorkData * workData,
                        int const vectorSz) const {
        auto * coeff = state.getReversedCoeffs().data();
        int const x_progress_earliest = x_and_ffts.progress - vectorSz + 1;
        for(int future=0; future<vectorSz; ++future) {
            int const x_unbounded_progress = x_progress_earliest + future;

            auto s = x_and_ffts.template getPastSegments<overlapMode>(x_unbounded_progress, n_coeffs);

            if(likely(s.size_from_start)) {
                T res;
                dotpr(&x_and_ffts.x[s.start],
                      coeff,
                      &res,
                      s.size_from_start);
                if constexpr (overlapMode == Overlap::Add) {
                    if(unlikely(s.size_from_zero)) {
                        T res2;
                        dotpr(&x_and_ffts.x[0],
                              coeff+s.size_from_start,
                              &res2,
                              s.size_from_zero);
                        res += res2;
                    }
                }
                y.writeOneFuture(future, res);
            }
        }
    }
    
    template<template<typename> typename Allocator2, typename WorkData>
    void step(State & state,
              XAndFFTS<T, Allocator2, Tag> const & x_and_ffts,
              Y<T, Tag> & y,
              WorkData * workData) const
    {
        auto s = x_and_ffts.template getPastSegments<overlapMode>(x_and_ffts.progress, n_coeffs);

        if(likely(s.size_from_start)) {
            using E = typename fft::RealSignal_<Tag, FPT>::type::value_type;
            T res;
            auto * coeff = state.getReversedCoeffs().data();
            dotpr(&x_and_ffts.x[s.start],
                  coeff,
                  &res,
                  s.size_from_start);
            if constexpr (overlapMode == Overlap::Add) {
                if(unlikely(s.size_from_zero)) {
                    E res2;
                    dotpr(&x_and_ffts.x[0],
                          coeff+s.size_from_start,
                          &res2,
                          s.size_from_zero);
                    res += res2;
                }
            }
            y.writeOneOutput(res);
        }
    }

    void flushToSilence(State & s) const {
    }
    
    int n_coeffs = 0;
};

}
