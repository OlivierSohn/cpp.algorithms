
namespace imajuscule {

struct AlgoScaleMetrics {
    int countCoeffs;
    int submissionPeriod;
};

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
    
    void logComputeState(Algo const & algo, std::ostream & os) const {
        os << "Custom scaling" << std::endl;
        IndentingOStreambuf indent(os);
        for(auto const & [param, algo] : v)
        {
            os << param.countCoeffs << " [" << param.submissionPeriod << "]" << std::endl;
            IndentingOStreambuf indent2(os);
            algo.logComputeState(os);
        }
    }
    
    bool isZero() const {
        return v.empty();
    }
    
    void flushToSilence() {
        for(auto & [param,c] : v) {
            c.flushToSilence();
        }
    }
    
    MinSizeRequirement setCoefficients(Algo const & algo,
                                       a64::vector<FPT> coeffs_)
    {
        v.clear();
        v.reserve(algo.v.size());
        
        auto it = coeffs_.begin();
        auto end = coeffs_.end();
        std::optional<MinSizeRequirement> res;

        int i=-1;
        for(auto & [param, conv] : algo.v) {
            ++i;
            if(it >= end) {
                // suboptimal CustomScaleConvolution
                break;
            }
            v.emplace_back();
            auto start = it;
            auto sizeBlock = param.countCoeffs;
            it += sizeBlock;
            MinSizeRequirement res2;
            if(it > end) {
                auto withPadding = a64::vector<FPT>{start,end};
                // We pad up-to sizeBlock. Benchmarks showed that this is time-wise better
                // than padding to the next power of 2 + delaying input.
                withPadding.resize(sizeBlock);
                res2 = v[i].setCoefficients(conv, std::move(withPadding));
            }
            else {
                res2 = v[i].setCoefficients(conv, {start,it});
            }
            
            if(!res) {
                res = res2;
            }
            else {
                res->mergeWith(res2);
            }
        }
        if(it < end) {
            throw std::runtime_error("not enough scales in CustomScaleConvolution");
        }

        if(res) {
            return *res;
        }
        return {0,0,0,{}};
    }

    double getEpsilon(Algo const & algo) const {
        double eps = 0.;
        
        Assert(v.size() <= algo.v.size());
        
        int i = -1;
        for(auto & state : v) {
            ++i;
            auto & [param, a] = algo.v[i];
            eps += state.getEpsilon(a);
        }
        
        return eps;
    }
    
    std::vector<A> v;
};

template<typename A>
struct AlgoCustomScaleConvolution {
    using SetupParam = CustomScaleConvolutionSetupParam<typename A::SetupParam>;

    using FPT = typename A::FPT;
    using Tag = typename A::Tag;
    using State = StateCustomScaleConvolution<typename A::State>;
    
    static constexpr auto zero_n_raw = fft::RealSignal_<Tag, FPT>::zero_n_raw;
    
    void setup(SetupParam const & p) {
        v.clear();
        v.reserve(p.scalingParams.size());
        
        for(auto const & param : p.scalingParams) {
            v.emplace_back();
            
            v.back().first.countCoeffs = param.countCoeffs;
            v.back().first.submissionPeriod = param.submissionPeriod;

            v.back().second.setup(param.setupParam);
            
            if(v.back().second.getLatency() != (param.submissionPeriod-1)) {
                // This breaks the class logic, and would lead to wrong results.
                throw std::logic_error("AlgoCustomScaleConvolution is applied to a type that doesn't respect the latency constraints.");
            }
        }
    }
    
    bool isValid() const {
        if(v.empty()) {
            return false;
        }
        return std::all_of(v.begin(), v.end(), [](auto & e) -> bool { return e.second.isValid(); });
    }

    void step(State & s,
              XAndFFTS<FPT, Tag> const & x_and_ffts,
              Y<FPT, Tag> & y)
    {
        Assert(s.v.size() <= v.size());
        
        int i = -1;
        for(auto & state : s.v) {
            ++i;
            auto & [param, algo] = v[i];

            // It is more optimal to do this check here than inside algo.doStep
            // because we early-exit at the first non matching parameter.

            auto const N = param.submissionPeriod;
            Assert(N > 0);
            Assert(is_power_of_two(N));
            
            if((x_and_ffts.fftsHalfSizesBitsToCompute & static_cast<uint32_t>(N)) != N) {
                // early exit
                return;
            }
            
            algo.forceStep(state,
                           x_and_ffts,
                           y);
        }
    }
    
    auto getLatency() const {
        if(v.empty()) {
            return 0;
        }
        return std::min_element(v.begin(),
                                v.end(),
                                [](auto const & p1, auto const & p2) {
            return p1.first.submissionPeriod < p2.first.submissionPeriod;
        })->first.submissionPeriod - 1;
    }
    
    int getWriteYBlockSize() const {
        return getBiggestScale();
    }
    
    int getBiggestScale() const {
        if(v.empty()) {
            return 0;
        }
        return std::max_element(v.begin(), v.end(), [](auto const & p1, auto const & p2){
            return p1.first.submissionPeriod < p2.first.submissionPeriod;
        })->first.submissionPeriod;
    }
    
public:
    std::vector<std::pair<AlgoScaleMetrics,A>> v;
};
}
