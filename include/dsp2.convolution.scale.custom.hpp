
namespace imajuscule {

template<typename A>
struct DescCustomScaleConvolution {
    static constexpr int nCoefficientsFadeIn = 0;
    static_assert(A::nCoefficientsFadeIn == 0); // else we need to handle them
    
    static constexpr bool has_subsampling = A::has_subsampling;
    static constexpr bool step_can_error = A::step_can_error;
    static_assert(!has_subsampling); // because it wouldn't make much sense
};

template<typename A>
struct AlgoCustomScaleConvolution;

template<typename A>
struct StateCustomScaleConvolution {
    using FPT = typename A::FPT;
    using Tag = typename A::Tag;
    using Desc = DescCustomScaleConvolution<typename A::Desc>;
    using Algo = AlgoCustomScaleConvolution<typename A::Algo>;
    
    static constexpr auto zero_n_raw = fft::RealSignal_<Tag, FPT>::zero_n_raw;
    
    void reset() {
        v.clear();
    }
    
    bool isZero() const {
        return v.empty();
    }
    
    void flushToSilence(Algo const & a) {
        auto itAlgo = a.v.begin();
        for(auto&c : v) {
            itAlgo->flushToSilence(c);
            ++itAlgo;
        }
    }

    static int getAllocationSz_SetCoefficients(typename Algo::SetupParam const & p) {
        int sz = 0;
        for(auto const & param : p.scalingParams) {
            sz += A::getAllocationSz_SetCoefficients(param.setupParam);
        }
        return sz;
    }
    void setCoefficients(Algo const & algo,
                         a64::vector<FPT> coeffs_)
    {
        v.clear();
        v.reserve(algo.v.size());
        
        auto it = coeffs_.begin();
        auto end = coeffs_.end();

        for(int i=0, endI=static_cast<int>(algo.v.size());
            i != endI;
            ++i)
        {
            if(it >= end) {
                // suboptimal
                break;
            }
            v.emplace_back();
            auto start = it;
            auto & conv = algo.v[i];
            if(i==endI-1) {
                it = end;
            }
            else {
                it += (algo.v[i+1].getLatency() - conv.getLatency()).toInteger();
            }
            v[i].setCoefficients(conv, {start,it});
        }
        if(it < end) {
            throw std::runtime_error("not enough scales in CustomScaleConvolution");
        }
    }
    
    template<typename F>
    void onContextFronteer(F f) {
        for(auto & e : v) {
            e.onContextFronteer(f);
        }
    }
    
    double getEpsilon(Algo const & algo) const {
        double eps = 0.;
        
        Assert(v.size() <= algo.v.size());
        
        int i = -1;
        for(auto & state : v) {
            ++i;
            eps += state.getEpsilon(algo.v[i]);
        }
        
        return eps;
    }

    void logComputeState(Algo const & algo, std::ostream & os) const {
        os << "Custom scaling" << std::endl;
        
        IndentingOStreambuf indent(os);
        
        for(int i=0, end=v.size(); i<end; ++i)
        {
            v[i].logComputeState(algo.v[i], os);
        }
    }
    
    std::vector<A> v;
};

template<typename A>
struct AlgoCustomScaleConvolution {
    template<typename TT>
    using Allocator = typename A::template Allocator<TT>;
    
    using SetupParam = CustomScaleConvolutionSetupParam<typename A::SetupParam>;
    using State = StateCustomScaleConvolution<typename A::State>;
    using Desc = DescCustomScaleConvolution<typename A::Desc>;

    using FPT = typename A::FPT;
    using Tag = typename A::Tag;
    
    static constexpr auto zero_n_raw = fft::RealSignal_<Tag, FPT>::zero_n_raw;
    
    
    static int getAllocationSz_Setup(SetupParam const & p) {
        int sz = 0;
        for(auto const & param : p.scalingParams) {
            sz += A::getAllocationSz_Setup(param.setupParam);
        }
        return sz;
    }
    void setup(SetupParam const & p) {
        v.clear();
        v.reserve(p.scalingParams.size());
        
        for(auto const & param : p.scalingParams) {
            v.emplace_back();
            v.back().setup(param.setupParam);
        }
    }
    
    bool handlesCoefficients() const {
        if(v.empty()) {
            return false;
        }
        return std::any_of(v.begin(), v.end(), [](auto & e) -> bool { return e.handlesCoefficients(); });
    }
    
    bool isValid() const {
        if(v.empty()) {
            return true;
        }
        return std::all_of(v.begin(), v.end(), [](auto & e) -> bool { return e.isValid(); });
    }

    void dephaseSteps(State & s,
                      int n_steps) const {
        // nothing. For justification, see comment in step()
    }

    template<template<typename> typename Allocator2>
    void step(State & s,
              XAndFFTS<FPT, Allocator2, Tag> const & x_and_ffts,
              Y<FPT, Tag> & y) const
    {
        Assert(s.v.size() <= v.size());
        
        int i = -1;
        
        // Note : the const is important here, and this is why the implementation of dephaseSteps is nop.
        for(auto const & state : s.v)
        {
            ++i;
            auto & algo = v[i];

            // It is more optimal to do this check here than inside algo.doStep
            // because we early-exit at the first non matching parameter.

            auto const N = algo.getBlockSize();
            Assert(N > 0);
            Assert(is_power_of_two(N));
            
            if((x_and_ffts.fftsHalfSizesBitsToCompute & static_cast<uint32_t>(N)) != N) {
                Assert(0 != (x_and_ffts.progress & static_cast<uint32_t>(N-1)));
                // early exit
                return;
            }
            // TODO no need to have fftsHalfSizesBitsToCompute if we can do this instead:
            Assert(0 == (x_and_ffts.progress & static_cast<uint32_t>(N-1)));
            
            algo.forceStep(state,
                           x_and_ffts,
                           y);
        }
    }
    
    void flushToSilence(State & s) const {
        s.flushToSilence(*this);
    }

    Latency getLatency() const {
        Assert(handlesCoefficients());
        return v.front().getLatency();
    }
    
    int getBiggestScale() const {
        if(v.empty()) {
            return 0;
        }
        return v.back().getBlockSize();
    }
    
public:
    std::vector<A> v;
};


template<typename A>
struct corresponding_legacy_dsp<AlgoCustomScaleConvolution<A>> {
    using type = CustomScaleConvolution<corresponding_legacy_dsp_t<A>>;
};

}
