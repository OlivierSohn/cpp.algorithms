
namespace imajuscule {

static constexpr int undefinedSplit = -1;
static constexpr int noSplit = -2;

/*
 * Creates a convolution scheme by combining 2 convolution schemes,
 * where A handles ** early ** coefficients, and B handles ** late ** coefficients.
 */
template<typename A, typename B>
struct SplitConvolution {
    using FPT = typename A::FPT;
    using EarlyHandler = A;
    using LateHandler = B;
    static constexpr int nComputePhaseable = A::nComputePhaseable + B::nComputePhaseable;
    static constexpr int nCoefficientsFadeIn = A::nCoefficientsFadeIn;
    static constexpr bool has_subsampling = LateHandler::has_subsampling; // we 'could' take earlyhandler into account, too.
    static constexpr bool step_can_error = A::step_can_error || B::step_can_error;

    struct SetupParam : public Cost {
        using AParam = typename EarlyHandler::SetupParam;
        using BParam = typename LateHandler::SetupParam;
        AParam aParams;
        BParam bParams;
        
        SetupParam(AParam a, BParam b)
        : aParams(a)
        , bParams(b)
        {}
        
        void logSubReport(std::ostream & os) const override {
            aParams.logSubReport(os);
            bParams.logSubReport(os);
        }
    };
    
    void logComputeState(std::ostream & os) const {
        
        os << "[";
        if(split==undefinedSplit) {
            os << "undefined";
        }
        else {
            os << "0..";
            if(split!=noSplit) {
                os << split-1;
            }
        }
        os << "]" << std::endl;
        
        {
            IndentingOStreambuf indent{os};
            a.logComputeState(os);
        }
        
        os << "[";
        if(split==undefinedSplit) {
            os << "undefined";
        } else {
            if(split!=noSplit) {
                os << split << "..";
            }
        }
        os << "]" << std::endl;
        
        {
            IndentingOStreambuf indent{os};
            b.logComputeState(os);
        }
    }
    
    static_assert(std::is_same_v<FPT,typename B::FPT>);
    
    bool isZero() const {
        return a.isZero() && b.isZero();
    }
    void reset() {
        a.reset();
        b.reset();
        split = undefinedSplit;
    }
    void flushToSilence() {
        a.flushToSilence();
        b.flushToSilence();
    }
    
    void setup(const SetupParam & p) {
        a.setup(p.aParams);
        b.setup(p.bParams);
        
        if(!b.isValid()) { // for case where lower resolution tail is not used
            split = noSplit;
        }
        else {
            split = B::nCoefficientsFadeIn + b.getLatency() - a.getLatency();
        }
    }
    
    void setCoefficients(a64::vector<FPT> coeffs_) {
        if(unlikely(undefinedSplit == split)) {
            throw std::logic_error("undefined split");
        }
        auto [rangeA,rangeB] = splitAt(split, coeffs_);
        
        if(rangeB.empty()) {
            a.setCoefficients(rangeA.materialize());
            b.reset();
            assert(b.isZero());
            return;
        }
        else {
            // prepend overlapping coefficients to b range:
            rangeB.setBegin(rangeB.begin() - B::nCoefficientsFadeIn);
            if(rangeB.begin() < rangeA.begin()) {
                LG(ERR, "sz_coeffs: %d, split: %d, nFadeIn : %d", coeffs_.size(), split, B::nCoefficientsFadeIn);
                throw std::logic_error("cannot cross-fade the coefficients");
            }
            a.setCoefficients(withLinearFadeOut(B::nCoefficientsFadeIn,rangeA.materialize()));
            b.setCoefficients(withLinearFadeIn( B::nCoefficientsFadeIn,rangeB.materialize()));
        }
        assert(!b.isZero());
        if(split + a.getLatency() == B::nCoefficientsFadeIn + b.getLatency()) {
            return;
        }
        LG(ERR, "split : %d, a lat: %d, b lat: %d, b init: %d", split, a.getLatency(), b.getLatency(), B::nCoefficientsFadeIn);
        throw std::logic_error("inconsistent split (latencies are not adapted)");
    }
    
    bool isValid() const {
        return a.isValid() && (b.isValid() || b.isZero());
    }
    
    template <typename Bool = bool>
    auto hasStepErrors() const -> std::enable_if_t<step_can_error, Bool> {
        bool errors = false;
        if constexpr (A::step_can_error) {
            errors = a.hasStepErrors();
        }
        if constexpr (B::step_can_error) {
            errors = b.hasStepErrors() || errors;
        }
        return errors;
    }
    FPT step(FPT val) {
        auto res = a.step(val);
        if(!b.isZero()) {
            res += b.step(val);
        }
        return res;
    }
    template<typename FPT2>
    void stepAddVectorized(FPT2 const * const input_buffer,
                           FPT2 * output_buffer,
                           int nSamples)
    {
        a.stepAddVectorized(input_buffer,
                            output_buffer,
                            nSamples);
        if(!b.isZero()) {
            b.stepAddVectorized(input_buffer,
                                output_buffer,
                                nSamples);
        }
    }
    template<typename FPT2>
    void stepAssignVectorized(FPT2 const * const input_buffer,
                              FPT2 * output_buffer,
                              int nSamples)
    {
        a.stepAssignVectorized(input_buffer,
                               output_buffer,
                               nSamples);
        if(!b.isZero()) {
            b.stepAddVectorized(input_buffer,
                                output_buffer,
                                nSamples);
        }
    }
    
    template<typename FPT2>
    void stepAddInputZeroVectorized(FPT2 * output_buffer,
                                    int nSamples)
    {
        a.stepAddInputZeroVectorized(output_buffer,
                                     nSamples);
        if(!b.isZero()) {
            b.stepAddInputZeroVectorized(output_buffer,
                                         nSamples);
        }
    }
    
    double getEpsilon() const {
        return a.getEpsilon() + b.getEpsilon() + std::numeric_limits<FPT>::epsilon();
    }
    
    auto getLatency() const {
        return a.getLatency();
    }
    
    std::array<int, nComputePhaseable> getComputePeriodicities() const {
        return arrayConcat(a.getComputePeriodicities(),
                           b.getComputePeriodicities());
    }
    // in [0, getComputePeriodicity())
    std::array<int, nComputePhaseable> getComputeProgresses() const {
        return arrayConcat(a.getComputeProgresses(),
                           b.getComputeProgresses());
    }
    void setComputeProgresses(std::array<int, nComputePhaseable> const & progresses) {
        auto [pa,pb] = arraySplit<A::nComputePhaseable,B::nComputePhaseable>(progresses);
        a.setComputeProgresses(pa);
        if(!b.isZero()) {
            b.setComputeProgresses(pb);
        }
    }
    
    auto & getA() { return a; }
    auto & getB() { return b; }
    auto const & getA() const { return a; }
    auto const & getB() const { return b; }
    
private:
    int split = undefinedSplit; // we have 'split' ** early ** coefficients, the rest are ** late ** coefficients.
    A a;
    B b;
};

template<typename A, typename B>
static std::ostream& operator << (std::ostream& os,
                                  const typename SplitConvolution<A,B>::SetupParam& p)
{
    using namespace std;
    os << "aParams:" << endl << p.aParams << endl << "bParams:" << endl << p.bParams << endl;
    return os;
}
    
}
