
namespace imajuscule {

template<typename A, typename B>
struct DescSplitConvolution {
    using LateHandler = B;

    static constexpr int nComputePhaseable = A::nComputePhaseable + B::nComputePhaseable;
    static constexpr int nCoefficientsFadeIn = A::nCoefficientsFadeIn;
    static constexpr bool has_subsampling = LateHandler::has_subsampling; // we 'could' take earlyhandler into account, too.
    static constexpr bool step_can_error = A::step_can_error || B::step_can_error;
};

template<typename A, typename B>
struct AlgoSplitConvolution;

template<typename A, typename B>
struct StateSplitConvolution {
    using Algo = AlgoSplitConvolution<typename A::Algo, typename B::Algo>;
    using Desc = DescSplitConvolution<typename A::Desc, typename B::Desc>;
    
    using FPT = typename A::FPT;
    using EarlyHandler = A;
    using LateHandler = B;
    
    static_assert(std::is_same_v<FPT,typename B::FPT>);
    
    bool isZero() const {
        return a.isZero() && b.isZero();
    }
    void reset() {
        a.reset();
        b.reset();
    }
    void flushToSilence() {
        a.flushToSilence();
        b.flushToSilence();
    }
    
    MinSizeRequirement setCoefficients(Algo const & algo, a64::vector<FPT> coeffs_) {
        int const split = algo.computeSplit();

        auto [rangeA,rangeB] = splitAt(split, coeffs_);

        if(rangeB.empty()) {
            auto res = a.setCoefficients(algo.getA(), rangeA.materialize());
            b.reset();
            assert(b.isZero());
            return res;
        }

        // prepend overlapping coefficients to b range:
        rangeB.setBegin(rangeB.begin() - B::Desc::nCoefficientsFadeIn);
        if(rangeB.begin() < rangeA.begin()) {
            LG(ERR, "sz_coeffs: %d, split: %d, nFadeIn : %d",
               coeffs_.size(),
               split,
               B::Desc::nCoefficientsFadeIn);
            throw std::logic_error("cannot cross-fade the coefficients");
        }
        auto res = a.setCoefficients(algo.getA(), withLinearFadeOut(B::Desc::nCoefficientsFadeIn,rangeA.materialize()));
        auto resB = b.setCoefficients(algo.getB(), withLinearFadeIn(B::Desc::nCoefficientsFadeIn,rangeB.materialize()));
        res.mergeWith(resB);
        assert(!b.isZero());
        return res;
    }
    
    template <typename Bool = bool>
    auto hasStepErrors() const -> std::enable_if_t<Desc::step_can_error, Bool> {
        bool errors = false;
        if constexpr (A::step_can_error) {
            errors = a.hasStepErrors();
        }
        if constexpr (B::step_can_error) {
            errors = b.hasStepErrors() || errors;
        }
        return errors;
    }
    
    double getEpsilon(Algo const & algo) const {
        return a.getEpsilon(algo.getA()) + b.getEpsilon(algo.getB()) + std::numeric_limits<FPT>::epsilon();
    }
    
    auto & getA() { return a; }
    auto & getB() { return b; }
    auto const & getA() const { return a; }
    auto const & getB() const { return b; }
    
private:
    A a;
    B b;
};

template<typename A, typename B>
struct AlgoSplitConvolution {
    using FPT = typename A::FPT;
    using EarlyHandler = A;
    using LateHandler = B;
    using State = StateSplitConvolution<typename A::State, typename B::State>;
    
    using Tag = typename A::Tag;
    static_assert(std::is_same_v<typename A::Tag, typename B::Tag>);

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
    
    static_assert(std::is_same_v<FPT,typename B::FPT>);
    
    void setup(const SetupParam & p) {
        a.setup(p.aParams);
        b.setup(p.bParams);
    }
    
    bool isValid() const {
        return a.isValid() && (b.isValid() /*|| b.isZero()*/); // il faudra peut-etre retravailler ca plus tard, pour finegrained
    }
    
    auto getLatency() const {
        return a.getLatency();
    }
    
    int computeSplit() const {
        if(!b.isValid()) { // for case where lower resolution tail is not used
            return noSplit;
        }
        else {
            return B::Desc::nCoefficientsFadeIn + b.getLatency() - a.getLatency();
        }
    }
    
    void step(State & s,
              XAndFFTS<FPT, Tag> const & x_and_ffts,
              Y<FPT, Tag> & y) const
    {
        a.step(s.getA(),
               x_and_ffts,
               y);
        b.step(s.getB(),
               x_and_ffts,
               y);
    }
    
    auto & getA() { return a; }
    auto & getB() { return b; }
    auto const & getA() const { return a; }
    auto const & getB() const { return b; }
    
private:
    A a;
    B b;
};

}
