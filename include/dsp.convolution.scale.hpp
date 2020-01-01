
namespace imajuscule {

struct CountDroppedScales {
    constexpr explicit CountDroppedScales(int n)
    : n(n)
    {}
    
    constexpr int toInteger() const {
        return n;
    }
private:
    int n;
};

struct ScaleConvolution_ {
    static constexpr int latencyForDroppedConvolutions(CountDroppedScales const & nDropped) {
        return static_cast<int>(pow2(static_cast<size_t>(nDropped.toInteger())))-1;
    }
    
    /*
     * The default value is optimal for the system I developped on (OSX / Macbook Air 2015 / intel core i7).
     * see 'TEST(BenchmarkConvolutions, scaled_vs_brute)' on how to determine the best value for a given system.
     *
     * TODO ideally this value should be a global, computed at initialization time.
     */
    static constexpr CountDroppedScales nDroppedOptimalFor_Split_Bruteforce_Fft = CountDroppedScales(6);
};

struct ScaleConvolutionSetupParam : public Cost {
    ScaleConvolutionSetupParam(CountDroppedScales nDropped = ScaleConvolution_::nDroppedOptimalFor_Split_Bruteforce_Fft)
    : nDropped(nDropped)
    {}
    
    /*
     * The number of early convolutions that are dropped.
     * The coefficients associated to these dropped convolutions are typically handled
     * by another convolution handler. (there are '2^nDropped - 1' such coefficients)
     *
     * (todo refactor to remove this warning : every PartitionAlgo should explicitely set this value)
     * WARNING: Do not change the default value unless you know what you are doing:
     * PartitionAlgo< ZeroLatencyScaledFineGrainedPartitionnedConvolution<T,FFTTag> >
     * assumes that this value will be the default value.
     */
    CountDroppedScales nDropped;
    
    int getImpliedLatency() const {
        return ScaleConvolution_::latencyForDroppedConvolutions(nDropped);
    }
    
    void logSubReport(std::ostream & os) const override {
        os << "Scaling, dropped: " << nDropped.toInteger() << std::endl;
    }
};

static std::vector<int> getScalingBlockSizes(int const szCoeffs, CountDroppedScales const & nDroppedConvolutions)
{
    assert( szCoeffs > 0 );
    
    std::vector<int> v;
    
    // note that these dropped coefficients are not passed to this function
    auto nFirstCoefficients = ScaleConvolution_::latencyForDroppedConvolutions(nDroppedConvolutions);
    int n = power_of_two_exponent(ceil_power_of_two(nFirstCoefficients+szCoeffs+1));
    Assert(nDroppedConvolutions.toInteger() < n);
    v.resize(n-nDroppedConvolutions.toInteger());
    for(int i=nDroppedConvolutions.toInteger(); i<n; ++i) {
        v[i-nDroppedConvolutions.toInteger()] = static_cast<int>(pow2(i));
    }
    return v;
}

template<typename A>
struct ScaleConvolutionSimulation {
    static constexpr bool has_subsampling = A::has_subsampling;
    static_assert(!has_subsampling); // because it wouldn't make much sense

    using FPT = typename A::FPT;
    using FFTTag = typename A::FFTTag;
    static constexpr int nComputePhaseable = 1;
    
    ScaleConvolutionSimulation()
    : costWriteOne(costWriteNConsecutive<FPT>(1))
    {}
    
    using SetupParam = ScaleConvolutionSetupParam;
    void setup(SetupParam const & p) {
        reset();
        nDroppedConvolutions = p.nDropped;
    }
    bool isZero() const {
        return v.empty();
    }
    void reset() {
        v.clear();
        x_halfSize = 0;
        nDroppedConvolutions = CountDroppedScales(0);
        reset_states();
    }
    
    void setCoefficientsCount(int64_t szCoeffs) {
        auto scalingSizes = getScalingBlockSizes(szCoeffs, nDroppedConvolutions);
        v.clear();
        reset_states();
        v.reserve(scalingSizes.size());
        
        for(auto sizeBlock : scalingSizes) {
            v.emplace_back();
            auto & conv = v.back();
            conv.setCoefficientsCount(sizeBlock);
            if(conv.getLatency() != sizeBlock - 1) {
                throw std::logic_error("Latency constraint violation");
            }
        }
        x_halfSize = scalingSizes.empty() ? 0 : scalingSizes.back();
    }
    
    auto getLatency() const {
        return ScaleConvolution_::latencyForDroppedConvolutions(nDroppedConvolutions);
    }
    std::array<int, nComputePhaseable> getComputePeriodicities() const {
        int res = x_halfSize;
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
        
        auto it = v.begin();
        auto end = v.end();
        
        double cost = costWriteOne;
        ++progress;
        endPadding = std::max(progress, endPadding);
        if(unlikely(x_halfSize == progress)) {
            Assert(endPadding == x_halfSize);
            // the second half of x is by design already zero padded
            endPadding = 2*x_halfSize;
        }
        
        int nUpdates = 1 + (static_cast<int>(count_trailing_zeroes(progress))) - nDroppedConvolutions.toInteger();
        if(nUpdates > 0) {
            auto endUpdate = it + nUpdates;
            
            for(int paddingSize = static_cast<int>(pow2(nDroppedConvolutions.toInteger()));
                it != endUpdate;
                ++it, paddingSize <<= 1)
            {
                // write padding
                int const neededEndPadding = progress + paddingSize;
                int const countPadding = neededEndPadding-endPadding;
                if(countPadding > 0) {
                    cost += fft::RealSignalCosts<FFTTag, FPT>::cost_zero_n_raw(countPadding);
                    endPadding = neededEndPadding;
                }
                cost += it->simuMajorStep();
            }
        }
        
        for(; it!= end; ++it) {
            cost += it->simuMinorStep();
        }
        
        if(unlikely(x_halfSize == progress)) {
            Assert(nUpdates == v.size());
            reset_states();
        }
        
        return cost;
    }
private:
    double costWriteOne;
    std::vector<A> v;
    int x_halfSize = 0, progress = 0, endPadding = 0;
    CountDroppedScales nDroppedConvolutions = CountDroppedScales(0);
    
    void reset_states() {
        progress = 0;
        endPadding = 0;
    }
};

/*
 * Creates a 0-latency convolution by combining 'n' sub-convolutions
 * of latencies 2^k-1, where k is in [0,n-1].
 *
 * The first convolutions can be dropped (see the constructor).
 */
template<typename A>
struct ScaleConvolution {
    using Simulation = ScaleConvolutionSimulation<typename A::Simulation>;
    using SetupParam = ScaleConvolutionSetupParam;
    
    using FPT = typename A::FPT;
    using RealSignal = typename A::RealSignal;
    using Tag = typename A::Tag;
    static constexpr auto zero_n_raw = fft::RealSignal_<Tag, FPT>::zero_n_raw;
    
    static constexpr int nComputePhaseable = 1;
    static constexpr int nCoefficientsFadeIn = 0;
    static_assert(A::nCoefficientsFadeIn == 0); // else we need to handle them
    
    static constexpr bool has_subsampling = A::has_subsampling;
    static constexpr bool step_can_error = A::step_can_error;
    static_assert(!has_subsampling); // because it wouldn't make much sense
    
    void logComputeState(std::ostream & os) const {
        os << "Scaling ["<< progress <<"/"<< x_halfSize <<"], dropped: " << nDroppedConvolutions.toInteger() << std::endl;
        IndentingOStreambuf indent(os);
        int i=nDroppedConvolutions.toInteger();
        for(auto const & algo : v)
        {
            ++i;
            os << i << std::endl;
            IndentingOStreambuf indent2(os);
            algo.logComputeState(os);
        }
    }
    
    bool isZero() const {
        return v.empty();
    }
    void setup(SetupParam const & p) {
        reset();
        nDroppedConvolutions = p.nDropped;
    }
    void reset() {
        v.clear();
        x.clear();
        x_halfSize = 0;
        nDroppedConvolutions = CountDroppedScales(0);
        reset_states();
    }
    void flushToSilence() {
        for(auto & c:v) {
            c.flushToSilence();
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
        auto scalingSizes = getScalingBlockSizes(coeffs_.size(), nDroppedConvolutions);
        v.clear();
        reset_states();
        v.reserve(scalingSizes.size());
        
        auto it = coeffs_.begin();
        auto end = coeffs_.end();
        for(auto sizeBlock : scalingSizes) {
            Assert(it <= end);
            auto start = it;
            it += sizeBlock;
            
            v.emplace_back();
            auto & conv = v.back();
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
            if(conv.getLatency() != sizeBlock - 1) {
                // This breaks the class logic, and would lead to wrong results.
                //   for example, FinegrainedPartitionnedFFTConvolution is not usable with this class.
                throw std::logic_error("Latency constraint violation");
            }
        }
        x_halfSize = scalingSizes.empty() ? 0 : scalingSizes.back();
        x.resize(2*x_halfSize); // including padding for biggest convolution

        progress = endPadding = 0;
        
        // fill the first half with a non-zero value to verify during tests that padding is done at the right time.
        std::fill(x.begin(),
                  x.begin() + x_halfSize,
                  typename RealSignal::value_type(1.9)
                  );
    }
    
    bool isValid() const {
        if(nDroppedConvolutions.toInteger() < 0) {
            return false;
        }
        if(v.empty()) {
            return true;
        }
        return std::all_of(v.begin(), v.end(), [](auto & e) -> bool { return e.isValid(); });
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
        return v[0].get_fft_length();
    }
private:
    FPT doStep(FPT val) {
        auto it = v.begin();
        auto end = v.end();
        
        x[progress] = typename RealSignal::value_type(val);
        ++progress;
        endPadding = std::max(progress, endPadding);
        if(unlikely(x_halfSize == progress)) {
            Assert(endPadding == x_halfSize);
            // the second half of x is by design already zero padded
            endPadding = 2*x_halfSize;
        }
        
        FPT r{};
        int nUpdates = 1 + (static_cast<int>(count_trailing_zeroes(progress))) - nDroppedConvolutions.toInteger();
        if(nUpdates > 0) {
            auto const xBeginPadding = x.begin() + progress;
            auto endUpdate = it + nUpdates;
            
            for(int paddingSize = static_cast<int>(pow2(nDroppedConvolutions.toInteger()));
                it != endUpdate;
                ++it, paddingSize <<= 1)
            {
                // write the padding exactly when we need it to optimize cache use
                {
                    int const neededEndPadding = progress + paddingSize;
                    int const countPadding = neededEndPadding - endPadding;
                    if(countPadding > 0) {
                        zero_n_raw(&x[endPadding], countPadding);
                        endPadding = neededEndPadding;
                    }
                }
                
                r += it->doStep(xBeginPadding - paddingSize);
            }
        }
        
        for(; it!= end; ++it) {
            r += it->doStep();
        }
        
        if(unlikely(x_halfSize == progress)) {
            Assert(nUpdates == v.size());
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
        return ScaleConvolution_::latencyForDroppedConvolutions(nDroppedConvolutions);
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
        return pow2(nDroppedConvolutions + v.size()) - pow2(nDroppedConvolutions);
    }
    
    int getBiggestScale() const {
        if(v.empty()) {
            return 0;
        }
        return pow2(nDroppedConvolutions.toInteger() + v.size() - 1);
    }
    
private:
    std::vector<A> v;
    RealSignal x;
    int x_halfSize = 0, progress = 0, endPadding = 0;
    CountDroppedScales nDroppedConvolutions = CountDroppedScales(0);

    void reset_states() {
        progress = 0;
        endPadding = 0;
    }
};

}
