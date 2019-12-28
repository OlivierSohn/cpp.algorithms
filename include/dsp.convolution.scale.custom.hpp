
namespace imajuscule {

template<typename SetupParam>
struct CustomScaleConvolutionSetupParam : public Cost {
    struct ScalingParam {
        int countCoeffs;
        int submissionPeriod;
        SetupParam setupParam;
    };

    CustomScaleConvolutionSetupParam(std::vector<ScalingParam> const & scalingParams)
    : scalingParams(scalingParams)
    {}
    
    std::vector<ScalingParam> scalingParams;
    
    int getImpliedLatency() const {
        if(scalingParams.empty()) {
            return 0;
        }
        return std::min_element(scalingParams.begin(),
                                scalingParams.end(),
                                [](auto const & p1, auto const & p2) {
            return p1.submissionPeriod < p2.submissionPeriod;
        })->submissionPeriod - 1;
    }
    
    void logSubReport(std::ostream & os) const override {
        os << "Custom scaling" << std::endl;
    }
};

struct ScaleMetrics {
    int countCoeffs;
    int submissionPeriod;
    int submissionCountdown;
};
struct CachedCosts {
    double minorCost;
    double majorCost;
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
            v.back().first.submissionPeriod = param.submissionPeriod;
            v.back().first.submissionCountdown = param.submissionPeriod-1;

            v.back().second.second.setup(param.setupParam);
        }
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
            param.submissionCountdown = param.submissionPeriod-1;
            auto sizeBlock = param.countCoeffs;
            algo.second.setCoefficientsCount(sizeBlock);
            totalCoeffs += sizeBlock;
            if(algo.second.getLatency() != param.submissionPeriod - 1) {
                // This breaks the class logic, and would lead to wrong results.
                throw std::logic_error("CustomScaleConvolutionSimulation is applied to a type that doesn't respect the latency constraints.");
            }
            // assuming these costs are constant
            algo.first.majorCost = algo.second.simuMajorStep();
            algo.first.minorCost = algo.second.simuMinorStep();
        }
        Assert(totalCoeffs < szCoeffs + v.back().first.countCoeffs); // account for possible padding
        x_halfSize = getBiggestScale();
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

    std::array<int, nComputePhaseable> getComputePeriodicities() const {
        int res = x_halfSize;
        Assert(res == getBiggestScale());
        return {res};
    }
    // in [0, getComputePeriodicity())
    std::array<int, nComputePhaseable> getComputeProgresses() const {
        return {static_cast<int>(progress)};
    }
    void setComputeProgresses(std::array<int, nComputePhaseable> const & progresses) {
        auto const p = progresses[0];
        while(getComputeProgresses()[0] != p) {
            simuStep();
        }
    }
    
    /*
     Dual method of ScaleConvolution::step()
     */
    double simuStep() {
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
                int paddingSize = param.submissionPeriod;
                // write the padding exactly when we need it to optimize cache use
                {
                    int const neededEndPadding = progress + paddingSize;
                    int const countPadding = neededEndPadding - endPadding;
                    if(countPadding > 0) {
                        cost += fft::RealSignalCosts<FFTTag, FPT>::cost_zero_n_raw(countPadding);
                        endPadding = neededEndPadding;
                    }
                }
                cost += algo.first.majorCost;
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
    std::vector<std::pair<ScaleMetrics,std::pair<CachedCosts, A>>> v;
    int x_halfSize = 0, progress = 0, endPadding = 0;
    
    void reset_states() {
        progress = 0;
        endPadding = 0;
    }
};


/*
 A generalization of ScaleConvolution
 */
template<typename A>
struct CustomScaleConvolution {
    using SetupParam = CustomScaleConvolutionSetupParam<typename A::SetupParam>;
    using Simulation = CustomScaleConvolutionSimulation<typename A::Simulation>;

    using FPT = typename A::FPT;
    using RealSignal = typename A::RealSignal;
    using Tag = typename A::Tag;
    static constexpr auto zero_n_raw = fft::RealSignal_<Tag, FPT>::zero_n_raw;
    
    static constexpr int nComputePhaseable = 1;
    static constexpr int nCoefficientsFadeIn = 0;
    static_assert(A::nCoefficientsFadeIn == 0); // else we need to handle them
    
    static constexpr bool has_subsampling = A::has_subsampling;
    static_assert(!has_subsampling); // because it wouldn't make much sense

    void logComputeState(std::ostream & os) const {
        os << "Custom scaling ["<< progress <<"/"<< x_halfSize <<"]" << std::endl;
        IndentingOStreambuf indent(os);
        for(auto const & [param, algo] : v)
        {
            os << param.countCoeffs << "[" << param.submissionCountDown << "/" << param.submissionPeriod << "]" << std::endl;
            IndentingOStreambuf indent2(os);
            algo.logComputeState(os);
        }
    }
    
    bool isZero() const {
        return v.empty();
    }
    void setup(SetupParam const & p) {
        reset();
        v.reserve(p.scalingParams.size());
        
        for(auto const & param : p.scalingParams) {
            v.emplace_back();
            
            v.back().first.countCoeffs = param.countCoeffs;
            v.back().first.submissionPeriod = param.submissionPeriod;
            v.back().first.submissionCountdown = param.submissionPeriod-1;

            v.back().second.setup(param.setupParam);
        }
    }
    void reset() {
        v.clear();
        x.clear();
        x_halfSize = 0;
        reset_states();
    }
    void flushToSilence() {
        for(auto & [param,c] : v) {
            c.flushToSilence();
            Assert(param.submissionPeriod > 0);
            param.submissionCountdown = param.submissionPeriod-1;
        }
        reset_states();
    }
    
    /*
     *   Unless the count of coefficients is of the form
     *     2^n - 1, some padding occurs.
     */
    // TODO we could require that 'A' tells what amount of underlying storage it will need, by number of coefficients,
    // and then this class could allocate the memory in one big chunk, and split it among 'A's.
    // This way, step / get could benefit from better memory locality, and memory accesses may be more predictable.
    void setCoefficients(a64::vector<FPT> coeffs_) {
        reset_states();
        
        auto it = coeffs_.begin();
        auto end = coeffs_.end();
        for(auto & [param, conv] : v) {
            param.submissionCountdown = param.submissionPeriod-1;
            Assert(it <= end);
            auto start = it;
            auto sizeBlock = param.countCoeffs;
            it += sizeBlock;
            if(it > end) {
                auto withPadding = a64::vector<FPT>{start,end};
                // We pad up-to sizeBlock. Benchmarks showed that this is time-wise better
                // than padding to the next power of 2 + delaying input.
                withPadding.resize(sizeBlock);
                conv.setCoefficients2(std::move(withPadding));
            }
            else {
                conv.setCoefficients2({start,it});
            }
            if(conv.getLatency() != param.submissionPeriod - 1) {
                // This breaks the class logic, and would lead to wrong results.
                throw std::logic_error("CustomScaleConvolution is applied to a type that doesn't respect the latency constraints.");
            }
        }
        Assert(it >= end);
        x_halfSize = getBiggestScale();
        x.resize(2*x_halfSize); // including padding for biggest convolution
        
        // fill the first half with a non-zero value to verify during tests that padding is done at the right time.
        std::fill(x.begin(),
                  x.begin() + x_halfSize,
                  typename RealSignal::value_type(1.9)
                  );
    }
    
    bool isValid() const {
        if(v.empty()) {
            return true;
        }
        return std::all_of(v.begin(), v.end(), [](auto & e) -> bool { return e.second.isValid(); });
    }
    
    /*
     Dual method of ScaleConvolutionSimulation::step()
     */
    FPT step(FPT val) {
        if(unlikely(isZero())) {
            return {};
        }
        return doStep(val);
    }
    
    int get_first_fft_length() const {
        if(isZero()) {
            return 0;
        }
        return v[0].second.get_fft_length();
    }
private:
    FPT doStep(FPT val) {
        x[progress] = typename RealSignal::value_type(val);
        ++progress;
        endPadding = std::max(progress, endPadding);
        if(unlikely(x_halfSize == progress)) {
            Assert(endPadding == x_halfSize);
            // the second half of x is by design already zero padded
            endPadding = 2*x_halfSize;
        }
        
        FPT r{};
        
        for(auto & [param, algo] : v) {
            Assert(param.submissionPeriod > 0);
            Assert(param.submissionCountdown >= 0);
            if(unlikely(0 == param.submissionCountdown)) {
                param.submissionCountdown = param.submissionPeriod;
                int paddingSize = param.submissionPeriod;
                // write the padding exactly when we need it to optimize cache use
                {
                    int const neededEndPadding = progress + paddingSize;
                    int const countPadding = neededEndPadding - endPadding;
                    if(countPadding > 0) {
                        zero_n_raw(&x[endPadding], countPadding);
                        endPadding = neededEndPadding;
                    }
                }
                r += algo.doStep(x.begin() + (progress - paddingSize));
            }
            else {
                r += algo.doStep();
            }
            --param.submissionCountdown;
        }
        
        if(unlikely(x_halfSize == progress)) {
            reset_states();
        }
        
        return r;
    }
public:
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
    void stepAssignVectorized(FPT2 const * const input_buffer,
                              FPT2 * output_buffer,
                              int nSamples)
    {
        if(unlikely(isZero())) {
            return;
        }
        for(int i=0; i<nSamples; ++i) {
            output_buffer[i] = doStep(input_buffer[i]);
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
        return epsilonOfNaiveSummation(v);
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

    std::array<int, nComputePhaseable> getComputePeriodicities() const {
        int res = x_halfSize;
        Assert(res == getBiggestScale());
        return {res};
    }
    // in [0, getComputePeriodicity())
    std::array<int, nComputePhaseable> getComputeProgresses() const {
        return {static_cast<int>(progress)};
    }
    void setComputeProgresses(std::array<int, nComputePhaseable> const & progresses) {
        auto const p = progresses[0];
        while(getComputeProgresses()[0] != p) {
            step(0);
        }
    }
    
    int countCoefficients() const {
        return std::accumulate(v.begin(), v.end(), 0, [](auto const & p) { return p.second.countCoefficients(); });
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
    std::vector<std::pair<ScaleMetrics,A>> v;
    RealSignal x;
    int x_halfSize = 0, progress = 0, endPadding = 0;
    
    void reset_states() {
        progress = 0;
        endPadding = 0;
    }
};

}
