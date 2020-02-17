
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
    
    static int getAllocationSz_SetCoefficients(typename Algo::SetupParam const & p) {
        return
        A::getAllocationSz_SetCoefficients(p.a) +
        B::getAllocationSz_SetCoefficients(p.b);
    }

    void setCoefficients(Algo const & algo, a64::vector<FPT> coeffs_) {
        int const split = algo.computeSplit();

        auto [rangeA,rangeB] = splitAt(split, coeffs_);

        if(rangeB.empty()) {
            a.setCoefficients(algo.a, rangeA.materialize());
            b.reset();
            return;
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
        a.setCoefficients(algo.a, withLinearFadeOut(B::Desc::nCoefficientsFadeIn,rangeA.materialize()));
        b.setCoefficients(algo.b, withLinearFadeIn(B::Desc::nCoefficientsFadeIn,rangeB.materialize()));
    }
    
    template<typename F>
    void onContextFronteer(F f) {
        a.onContextFronteer(f);
        b.onContextFronteer(f);
    }
    
    template <typename Bool = bool>
    auto hasStepErrors() const -> std::enable_if_t<Desc::step_can_error, Bool> {
        bool errors = false;
        if constexpr (A::Desc::step_can_error) {
            errors = a.hasStepErrors();
        }
        if constexpr (B::Desc::step_can_error) {
            errors = b.hasStepErrors() || errors;
        }
        return errors;
    }
    
    double getEpsilon(Algo const & algo) const {
        return a.getEpsilon(algo.a) + b.getEpsilon(algo.b) + std::numeric_limits<FPT>::epsilon();
    }
    void logComputeState(Algo const & algo, std::ostream & os) const {
        os << "[";
        int const split = algo.computeSplit();
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
            a.logComputeState(algo.a, os);
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
            b.logComputeState(algo.b, os);
        }
    }

    A a;
    B b;
};

template<typename A, typename B>
struct AlgoSplitConvolution {
    template<typename TT>
    using Allocator = typename A::template Allocator<TT>;

    // see https://stackoverflow.com/questions/38630445/stdis-same-equivalent-for-unspecialised-template-types
    /*
    template<typename TT>
    using AllocatorB = typename B::template Allocator<TT>;
    static_assert(std::is_same_v<Allocator, AllocatorB>);
     */

    using FPT = typename A::FPT;
    using EarlyHandler = A;
    using LateHandler = B;
    using State = StateSplitConvolution<typename A::State, typename B::State>;
    using Desc = DescSplitConvolution<typename A::Desc, typename B::Desc>;
    
    using Tag = typename A::Tag;
    static_assert(std::is_same_v<typename A::Tag, typename B::Tag>);

    using SetupParam = SplitSetupParam<typename A::SetupParam, typename B::SetupParam>;

    static_assert(std::is_same_v<FPT,typename B::FPT>);
    
    void setup(const SetupParam & p) {
        a.setup(p.a);
        b.setup(p.b);
    }

    bool isValid() const {
        return a.isValid() && b.isValid();
    }
    bool handlesCoefficients() const {
        return a.handlesCoefficients();
    }

    Latency getLatency() const {
        Assert(handlesCoefficients());
        return a.getLatency();
    }
    
    int getBiggestScale() const {
        return std::max(a.getBiggestScale(),
                        b.getBiggestScale());
    }
    
    int computeSplit() const {
        if(!b.handlesCoefficients()) {
            return noSplit;
        }
        else {
            return B::Desc::nCoefficientsFadeIn + (b.getLatency() - a.getLatency()).toInteger();
        }
    }

    void dephaseStep(State &s,
                     int x_progress) const {
        a.dephaseStep(s.a, x_progress);
        b.dephaseStep(s.b, x_progress);
    }
    
    template<template<typename> typename Allocator2, typename WorkData>
    void stepVectorized(State & s,
                        XAndFFTS<FPT, Allocator2, Tag> const & x_and_ffts,
                        Y<FPT, Tag> & y,
                        WorkData * workData,
                        int const vectorSz) const {
        // We step the latehandler before the early handler to minimize
        // the number of splits done in Y::addAssign and Y::addAssignPresent
        // (latehandlers write bigger chunks than early handlers)
        b.stepVectorized(s.b,
                         x_and_ffts,
                         y,
                         workData,
                         vectorSz);
        a.stepVectorized(s.a,
                         x_and_ffts,
                         y,
                         workData,
                         vectorSz);
    }
    
    template<template<typename> typename Allocator2, typename WorkData>
    void step(State & s,
              XAndFFTS<FPT, Allocator2, Tag> const & x_and_ffts,
              Y<FPT, Tag> & y,
              WorkData * workData) const
    {
        // We step the latehandler before the early handler to minimize
        // the number of splits done in Y::addAssign and Y::addAssignPresent
        // (latehandlers write bigger chunks than early handlers)
        b.step(s.b,
               x_and_ffts,
               y,
               workData);
        a.step(s.a,
               x_and_ffts,
               y,
               workData);
    }
    
    void flushToSilence(State & s) const {
        a.flushToSilence(s.a);
        b.flushToSilence(s.b);
    }

    A a;
    B b;
};


}
