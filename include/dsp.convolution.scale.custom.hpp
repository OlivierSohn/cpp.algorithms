
namespace imajuscule {

template<typename SetupParam>
struct CustomScaleConvolutionSetupParam : public Cost {
    
    static constexpr bool has_subsampling = false;
    static_assert(!SetupParam::has_subsampling);
    static constexpr int nCoefficientsFadeIn = 0;
    static_assert(SetupParam::nCoefficientsFadeIn == 0); // else we need to handle them

    struct ScalingParam {
        ScalingParam(int countCoeffs,
                     SetupParam setupParam)
        : countCoeffs(countCoeffs)
        , setupParam(setupParam)
        {}
        
        int countCoeffs;
        SetupParam setupParam;
        
        bool isValid() const {
            if(countCoeffs < 0) {
                return false;
            }
            return setupParam.isValid();
        }
        bool handlesCoefficients() const {
            return countCoeffs > 0;
        }
        
        void logSubReport(std::ostream & os) const {
            os << countCoeffs << " coeffs with:" << std::endl;

            {
                IndentingOStreambuf indent(os);
                setupParam.logSubReport(os);
            }
        }
    };

    CustomScaleConvolutionSetupParam(std::vector<ScalingParam> const & scalingParams)
    : scalingParams(scalingParams)
    {}
    
    std::vector<ScalingParam> scalingParams;
    
    void adjustWork(int const nTargetCoeffs) {
        int nHandledCoeffs = 0;
        int sz = 0;
        for(auto & s : scalingParams) {
            if(nHandledCoeffs >= nTargetCoeffs) {
                break;
            }
            int const nRemainingCoeffs = nTargetCoeffs - nHandledCoeffs;
            Assert(s.countCoeffs <= s.setupParam.countMaxHandledCoeffs());
            if(nRemainingCoeffs < s.setupParam.countMaxHandledCoeffs()) {
                s.setupParam.adjustWork(nRemainingCoeffs);
            }
            nHandledCoeffs += s.setupParam.countMaxHandledCoeffs();
            ++sz;
        }
        
        Assert(nHandledCoeffs >= nTargetCoeffs); // else, we need to expand the last scale

        // alternative to std::vector::resize (because ScalingParam has no default constructor)
        while(scalingParams.size() > sz) {
            scalingParams.pop_back();
        }
        if(!handlesCoefficients()) {
            setCost(0.);
        }
    }
    
    template<Overlap Mode, typename FFTAlgo>
    MinSizeRequirement getMinSizeRequirement(int const maxVectorSz) const
    {
        std::optional<MinSizeRequirement> res;
        for(auto const & a : scalingParams)
        {
            auto res2 = a.setupParam.template getMinSizeRequirement<Mode, FFTAlgo>(maxVectorSz);
            if(!res) {
                res = res2;
            }
            else {
                res->mergeWith(res2);
            }
        }
        if(res) {
            return *res;
        }
        return {0,0,{},0};
    }

    Latency getImpliedLatency() const {
        Assert(handlesCoefficients());
        return scalingParams.begin()->setupParam.getImpliedLatency();
    }
    
    int getBiggestScale() const {
        int res = 0;
        for(auto const & p : scalingParams) {
            res = std::max(res,
                           p.setupParam.getBiggestScale());
        }
        return res;
    }
    
    void logSubReport(std::ostream & os) const override {
        os << "Custom scaling" << std::endl;

        IndentingOStreambuf indent(os);
        for(auto const & s : scalingParams) {
            s.logSubReport(os);
        }
    }
    
    bool isValid() const {
        if(scalingParams.empty()) {
            return true;
        }
        return std::all_of(scalingParams.begin(),
                           scalingParams.end(),
                           [](auto const & p) -> bool { return p.isValid(); });
    }
    bool handlesCoefficients() const {
        if(scalingParams.empty()) {
            return false;
        }
        return scalingParams.begin()->handlesCoefficients();
    }

    template<typename F>
    void forEachUsingSameContext(F f) const {
        f(*this);
        for(auto & s : scalingParams) {
            s.setupParam.forEachUsingSameContext(f);
        }
    }
};


struct SimuScaleMetrics {
    int countCoeffs;
    int submissionPeriod;
    int submissionCountdown;
};
struct ScaleMetrics {
    int countCoeffs;
    int submissionCountdown;
};
struct CachedCosts {
    double minorCost;
};

template<typename A>
struct CustomScaleConvolutionSimulation {
    static constexpr bool has_subsampling = A::has_subsampling;
    static_assert(!has_subsampling); // because it wouldn't make much sense

    using FPT = typename A::FPT;
    using FFTTag = typename A::FFTTag;
    static constexpr int nComputePhaseable = 1;
    
    using SetupParam = CustomScaleConvolutionSetupParam<typename A::SetupParam>;
    
    CustomScaleConvolutionSimulation()
    : costWriteOne(costWriteNConsecutive<FPT>(1))
    {}
    
    bool isZero() const {
        return v.empty();
    }
    void setup(SetupParam const & p) {
        reset();
        v.reserve(p.scalingParams.size());
        
        for(auto const & param : p.scalingParams) {
            v.emplace_back();
            
            v.back().first.countCoeffs = param.countCoeffs;
            v.back().first.submissionPeriod = param.setupParam.getImpliedLatency().toInteger() + 1;
            v.back().first.submissionCountdown = v.back().first.submissionPeriod - 1;

            v.back().second.second.setup(param.setupParam);
        }

        x_halfSize = getBiggestScale();
    }
    void reset() {
        v.clear();
        x_halfSize = 0;
        reset_states();
    }
    
    void setCoefficientsCount(int64_t szCoeffs) {
        reset_states();
        
        int totalCoeffs = 0;
        for(auto & [param, algo] : v) {
            auto sizeBlock = param.countCoeffs;
            algo.second.setCoefficientsCount(sizeBlock);
            totalCoeffs += sizeBlock;
            param.submissionCountdown = param.submissionPeriod-1;
            if(algo.second.getLatency().toInteger() != param.submissionCountdown) {
                throw std::logic_error("CustomScaleConvolutionSimulation is applied to a type that doesn't respect the latency constraints.");
            }
            // assuming these costs are constant
            algo.first.minorCost = algo.second.simuMinorStep();
        }
        if(szCoeffs > totalCoeffs) {
            throw std::logic_error("too much coefficients for scales");
        }
        if(!v.empty()) {
            if(szCoeffs < totalCoeffs - v.back().first.countCoeffs + 1) {
                throw std::logic_error("not enough coefficients for scales");
            }
        }
    }

    void dephaseByGroupRatio(float phase_group_ratio)
    {
        Assert(phase_group_ratio >= 0.f);
        Assert(phase_group_ratio < 1.f);

        int const phase_period = x_halfSize;
        Assert(phase_period == getBiggestScale());
        
        XFFTsCostsFactors dummy; // not important because we don't use the return value of simuStep.
        int nsteps = static_cast<int>(0.5f + phase_group_ratio * phase_period);
        for(int i=0; i<nsteps; ++i) {
            simuStep(dummy);
        }
    }
    
    /*
     Dual method of CustomScaleConvolution::step()
     */
    double simuStep(XFFTsCostsFactors const & xFftCostFactors) {
        if(unlikely(isZero())) {
            return {};
        }
        
        double cost = costWriteOne;
        ++progress;
        endPadding = std::max(progress, endPadding);
        if(unlikely(x_halfSize == progress)) {
            Assert(endPadding == x_halfSize);
            // the second half of x is by design already zero padded
            endPadding = 2*x_halfSize;
        }
                
        for(auto & [param, algo] : v) {
            Assert(param.submissionPeriod > 0);
            Assert(param.submissionCountdown >= 0);
            if(unlikely(0 == param.submissionCountdown)) {
                param.submissionCountdown = param.submissionPeriod;
                int const paddingSize = param.submissionPeriod;
                // write the padding exactly when we need it to optimize cache use
                {
                    int const neededEndPadding = progress + paddingSize;
                    int const countPadding = neededEndPadding - endPadding;
                    if(countPadding > 0) {
                        cost += fft::RealSignalCosts<FFTTag, FPT>::cost_zero_n_raw(countPadding);
                        endPadding = neededEndPadding;
                    }
                }
                cost += algo.second.simuMajorStep(xFftCostFactors);
            }
            else {
                cost += algo.first.minorCost;
            }
            --param.submissionCountdown;
        }
        
        if(unlikely(x_halfSize == progress)) {
            reset_states();
        }
        
        return cost;
    }
    
    double simuBatch(int64_t nRemainingSteps,
                     XFFTsCostsFactors const & xFftCostFactors) {
        double cost{};
        
        while(progress != 0 && nRemainingSteps) {
            //LG(INFO, "! pre");
            Assert(0); // by design a batch should start at progress == 0
            cost += simuStep(xFftCostFactors);
            --nRemainingSteps;
        }
        
        int64_t const period = getBiggestScale();
        //LG(INFO, "period %d", period);
        if(period) {
            int64_t const nFullPeriods = nRemainingSteps / period;
            nRemainingSteps -= nFullPeriods * period;
            //LG(INFO, "%d periods", nFullPeriods);
            cost += nFullPeriods * simuBiggestPeriod(xFftCostFactors);
        }
        
        while(nRemainingSteps) {
            //LG(INFO, "! post");
            Assert(0); // by design a batch should be a full number of periods
            cost += simuStep(xFftCostFactors);
            --nRemainingSteps;
        }
        //LG(INFO, "%f", cost);
        return cost;
    }

    int getBiggestScale() const {
        if(v.empty()) {
            return 0;
        }
        return std::max_element(v.begin(), v.end(), [](auto const & p1, auto const & p2){
            return p1.first.submissionPeriod < p2.first.submissionPeriod;
        })->first.submissionPeriod;
    }
private:
    double costWriteOne;
    std::vector<std::pair<SimuScaleMetrics,std::pair<CachedCosts, A>>> v;
    int x_halfSize = 0, progress = 0, endPadding = 0;
    
    void reset_states() {
        progress = 0;
        endPadding = 0;
    }
    
    double simuBiggestPeriod(XFFTsCostsFactors const & xFftCostFactors) {
        // we compute the cost for 'biggestSubmissionPeriod' iterations
        int const biggestSubmissionPeriod = getBiggestScale();
        
        double cost = costWriteOne * biggestSubmissionPeriod;

        // padding cost
        if(v.size() >= 2)
        {
            int Gcd = gcd(v[0].first.submissionPeriod,
                          v[1].first.submissionPeriod);
            for(auto it = v.begin() + 2,
                end = v.end();
                it < end;
                ++it)
            {
                Gcd = gcd(Gcd, it->first.submissionPeriod);
            }
            Assert(Gcd > 0);
            
            int endPadding = 0;
            for(int progress = Gcd;
                progress < biggestSubmissionPeriod;
                progress += Gcd)
            {
                endPadding = std::max(progress, endPadding);

                for(auto & [param, algo] : v) {
                    Assert(param.submissionPeriod > 0);
                    if(param.submissionPeriod == biggestSubmissionPeriod) {
                        // no padding for this one
                        continue;
                    }
                    if(progress * (param.submissionPeriod / progress) == param.submissionPeriod) {
                        int const paddingSize = param.submissionPeriod;
                        int const neededEndPadding = progress + paddingSize;
                        int const countPadding = neededEndPadding - endPadding;
                        if(countPadding > 0) {
                            cost += fft::RealSignalCosts<FFTTag, FPT>::cost_zero_n_raw(countPadding);
                            endPadding = neededEndPadding;
                        }
                    }
                }
            }
        }
        
        //LG(INFO, "  padding %f", cost);
        
        // major / minor costs
        for(auto & [param, algo] : v) {
            Assert(param.submissionPeriod > 0);
            int nMajorSteps = biggestSubmissionPeriod / param.submissionPeriod;
            Assert(nMajorSteps * param.submissionPeriod == biggestSubmissionPeriod);

            int nMinorSteps = biggestSubmissionPeriod - nMajorSteps;
            Assert(nMinorSteps >= 0);
            
            //LG(INFO, "  param.submissionPeriod %d", param.submissionPeriod);
            //LG(INFO, "  %d major steps", nMajorSteps);
            cost += nMajorSteps * algo.second.simuMajorStep(xFftCostFactors);
            //LG(INFO, "  %d minor steps", nMinorSteps);
            cost += nMinorSteps * algo.first.minorCost;
        }
        return cost;
    }
};

template<typename A, typename T, typename FFTTag>
struct Simulation_<CustomScaleConvolutionSetupParam<A>, T, FFTTag> {
    using type = CustomScaleConvolutionSimulation<Simulation<A, T, FFTTag>>;
};


}
