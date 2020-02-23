
namespace imajuscule {

static constexpr int undefinedSplit = -1;
static constexpr int noSplit = -2;

template<typename A, typename B>
struct SplitSetupParam : public Cost {
    static constexpr bool has_subsampling = B::has_subsampling;
    static constexpr int nCoefficientsFadeIn = A::nCoefficientsFadeIn;

    using AParam = A;
    using BParam = B;
    AParam a;
    BParam b;
    
    SplitSetupParam(AParam const & a,
                    BParam const & b)
    : a(a)
    , b(b)
    {}
    
    void logSubReport(std::ostream & os) const override {
        a.logSubReport(os);
        b.logSubReport(os);
    }
    
    bool handlesCoefficients() const {
        return a.handlesCoefficients();
    }
    
    Latency getImpliedLatency() const {
        Assert(handlesCoefficients());
        return a.getImpliedLatency();
    }

    int getBiggestScale() const {
        return std::max(a.getBiggestScale(),
                        b.getBiggestScale());
    }

    template<typename F>
    void forEachUsingSameContext(F f) const {
        f(*this);
        a.forEachUsingSameContext(f);
        b.forEachUsingSameContext(f);
    }
    
    void adjustWork(int targetNCoeffs) {
        if(!b.handlesCoefficients()) {
            a.adjustWork(targetNCoeffs);
        }
        else {
            int const split = B::nCoefficientsFadeIn + (b.getImpliedLatency() - a.getImpliedLatency()).toInteger();
            Assert(split >= 0);
            int const nEarly = std::min(split, targetNCoeffs);
            a.adjustWork(nEarly);
            b.adjustWork(std::max(0,
                                  targetNCoeffs - nEarly + B::nCoefficientsFadeIn));
        }
        if(!handlesCoefficients()) {
            setCost(0.);
        }
    }
    
    template<Overlap Mode, typename FFTAlgo>
    MinSizeRequirement getMinSizeRequirement(int const maxVectorSz) const
    {
        auto res = a.template getMinSizeRequirement<Mode, FFTAlgo>(maxVectorSz);
        res.mergeWith(b.template getMinSizeRequirement<Mode, FFTAlgo>(maxVectorSz));
        return res;
    }
};


template<typename A, typename B, typename T, typename Tag>
struct SplitConvolutionSimulation {
    using SimA = Simulation<A, T, Tag>;
    using SimB = Simulation<B, T, Tag>;
    
    void setup(SplitSetupParam<A,B> const & p) {
        a.setup(p.a);
        b.setup(p.b);
    }
    double simuStep(XFFTsCostsFactors const & xFFTCostFactors) {
        return b.simuStep(xFFTCostFactors) + a.simuStep(xFFTCostFactors);
    }
    
    int getBiggestScale() const {
        int const sa = a.getBiggestScale();
        int const sb = b.getBiggestScale();
        if(!sb) {
            return sa;
        }
        if(!sa) {
            return sb;
        }
        int64_t const p = ppcm(sa, sb);
        if(p > std::numeric_limits<int>::max()) {
            throw std::runtime_error("biggest scale overflows int");
        }
        Assert(p == std::max(sa, sb)); // out of curiosity
        return static_cast<int>(p);
    }
private:
    SimA a;
    SimB b;
};

template<typename A, typename B, typename T, typename Tag>
struct Simulation_<SplitSetupParam<A, B>, T, Tag> {
    using type = SplitConvolutionSimulation<A, B, T, Tag>;
};

template<typename A, typename B>
static std::ostream& operator << (std::ostream& os,
                                  const SplitSetupParam<A,B>& p)
{
    using namespace std;
    os << "a:" << endl << p.a << endl << "b:" << endl << p.b << endl;
    return os;
}
    
}
