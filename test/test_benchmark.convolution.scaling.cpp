namespace imajuscule {


template<typename T, template<typename> typename Allocator, typename Tag>
auto mkConvolution(std::vector<Scaling> const & v,
                   a64::vector<T> const & coeffs) {
    using CLegacy = CustomScaleConvolution<FFTConvolutionIntermediate < PartitionnedFFTConvolutionCRTP<T, Allocator, Tag> >>;
    using CNew = Convolution<
    AlgoCustomScaleConvolution<AlgoFFTConvolutionIntermediate < AlgoPartitionnedFFTConvolutionCRTP<T, Allocator, Tag> >>
    >;
    //using C = CLegacy;
    using C = CNew;
    
    using SetupParam = typename C::SetupParam;
    using ScalingParam = typename SetupParam::ScalingParam;
    
    auto scalingParams = scalingsToParams<ScalingParam>(v);
    auto c = std::make_unique<C>();
    
    SetupParam p({scalingParams});

    int const sz =
    C::getAllocationSz_Setup(p) +
    C::getAllocationSz_SetCoefficients(p);
    
    using FBinsAllocator = typename fft::RealFBins_<Tag, T, Allocator>::type::allocator_type;
    
    using MemRsc = MemResource<FBinsAllocator>;
    
    auto memory = std::make_unique<typename MemRsc::type>(sz);
    
    {
        typename MemRsc::use_type use(*memory);

        c->setup(p);
        c->setCoefficients(coeffs);
    }
    
    if constexpr (MemRsc::limited) {
        //std::cout << "mem monotonic " << memory->used() << std::endl;
        Assert(0 == memory->remaining());
        // else we should fix getAllocationSz_Setup / getAllocationSz_SetCoefficients
    }

    return std::make_pair(std::move(c), std::move(memory));
}


struct Costs
{
    double simulated, real, real_flushes;
    
    Costs operator * (double v) const {
        return {simulated*v, real*v, real_flushes*v};
    }
};

template<typename T, template<typename> typename Allocator, typename Tag>
void findCheapest2(int const firstSz, int const nCoeffs, XFFtCostFactors const & factors, double & sideEffect)
{
    using C = CustomScaleConvolution<FFTConvolutionIntermediate < PartitionnedFFTConvolutionCRTP<T, Allocator, Tag> >>;
    
    auto coeffs = mkTestCoeffs<T>(nCoeffs);
    
    std::map<std::vector<Scaling>, Costs> results;
    {
        int count = 0;
        ScalingsIterator{
            firstSz,
            static_cast<int>(coeffs.size())
        }.forEachScaling([&coeffs, &count, &results, &factors, &sideEffect](auto const & v){
            auto [conv, mem] = mkConvolution<T, Allocator, Tag>(v, coeffs);
            constexpr int cache_flush_period_samples = 128;
            auto sim = mkSimulation<C>(v, coeffs.size());
            results.emplace(v,
                            Costs{
                virtualCostPerSample(sim, factors),
                realCostPerSample(*conv, sideEffect),
                realCostPerSampleWithCacheFlushes(*conv, cache_flush_period_samples, sideEffect)
            });
            ++count;
        });
        if(results.size() != count) {
            throw std::logic_error("duplicate results");
        }
    }
    auto bySimulatedCost = orderByValue(results, [](auto it1, auto it2){
        return it1->second.simulated < it2->second.simulated;
    });
    auto byRealCost = orderByValue(results, [](auto it1, auto it2){
        return it1->second.real < it2->second.real;
    });
    auto byRealFlushesCost = orderByValue(results, [](auto it1, auto it2){
        return it1->second.real_flushes < it2->second.real_flushes;
    });
    auto & bestByRealCost = byRealCost[0];
    auto & bestByRealFlushesCost = byRealFlushesCost[0];
    auto & bestBySimulatedCost = bySimulatedCost[0];
    if(bestByRealCost != bestBySimulatedCost || bestByRealCost != bestByRealFlushesCost) {
        std::cout << "mismatch" << std::endl;
        displayScaling(bestByRealFlushesCost, firstSz, nCoeffs);
        displayScaling(bestByRealCost, firstSz, nCoeffs);
        displayScaling(bestBySimulatedCost, firstSz, nCoeffs);
    }
    else {
        std::cout << "==" << std::endl;
        displayScaling(bestByRealCost, firstSz, nCoeffs);
    }
    
    /*
     int const nDisplay = -1;
     analyzeScalings(firstSz, nCoeffs, bySimulatedCost, nDisplay);
     */
}
template<typename T, template<typename> typename Allocator, typename Tag>
void findCheapest2()
{
    XFFtCostFactors factors;
    double sideEffect{};
    for(int factor = 100; factor <= 10000; factor *=10)
    {
        for(int firstSz = 1; firstSz < 100000; firstSz*=2)
        {
            int const nCoeffs = factor*firstSz;
            if(nCoeffs > 2000000) {
                continue;
            }
            findCheapest2<T, Allocator, Tag>(firstSz, nCoeffs, factors, sideEffect);
        }
    }
    std::cout << "sideEffect " << sideEffect << std::endl;
}


template<typename T, template<typename> typename Allocator>
void findCheapest()
{
    for_each(fft::Tags, [](auto t) {
        using Tag = decltype(t);
        COUT_TYPE(Tag); std::cout << std::endl;
        findCheapest2<T, Allocator, Tag>();
    });
}
enum class CostModel {
    VirtualSimulation,
    RealSimulation,
    RealSimulationCacheFlushes
};

template<typename T, template<typename> typename Allocator, typename Tag>
void analyzeSimulated(int const firstSz, int const nCoeffs, CostModel model, XFFtCostFactors const & factors, double & sideEffect)
{
    using C = CustomScaleConvolution<FFTConvolutionIntermediate < PartitionnedFFTConvolutionCRTP<T, Allocator, Tag> >>;
    
    auto coeffs = mkTestCoeffs<T>(nCoeffs);
    
    std::map<std::vector<Scaling>, double> results;
    {
        int count = 0;
        ScalingsIterator{
            firstSz,
            static_cast<int>(coeffs.size())
        }.forEachScaling([model, &coeffs, &count, &results, &factors, &sideEffect](auto const & v){
            //displayScalingNumerically(v, std::cout); std::cout << std::endl;
            /*if(v.size() != 1) {
                return;
            }*/
            switch(model) {
                case CostModel::VirtualSimulation:
                {
                    auto sim = mkSimulation<C>(v, coeffs.size());
                    results.emplace(v,
                                    virtualCostPerSample(sim, factors));
                }
                    break;
                case CostModel::RealSimulation:
                {
                    auto [conv, mem] = mkConvolution<T, Allocator, Tag>(v, coeffs);
                    results.emplace(v,
                                    realCostPerSample(*conv, sideEffect));
                    
                }
                    break;
                case CostModel::RealSimulationCacheFlushes:
                {
                    auto [conv, mem] = mkConvolution<T, Allocator, Tag>(v, coeffs);
                    constexpr int cache_flush_period_samples = 128;
                    results.emplace(v,
                                    realCostPerSampleWithCacheFlushes(*conv, cache_flush_period_samples, sideEffect));
                }
                    break;
            }
            ++count;
        });
        if(results.size() != count) {
            throw std::logic_error("duplicate results");
        }
    }
    auto bySimulatedCost = orderByValue(results, [](auto it1, auto it2){
        return it1->second < it2->second;
    });
    //            auto & bestBySimulatedCost = bySimulatedCost[0];
    //            displayScaling(bestBySimulatedCost, firstSz, nCoeffs);
    
    int const nDisplay = -1;
    analyzeScalings(firstSz, nCoeffs, bySimulatedCost, nDisplay);
}

template<typename T, template<typename> typename Allocator>
void smallTest(void) {
    double sideEffect{};
    int const scale = 1;
    for(int i=0; i<3; ++i) {
        CostModel costmodel;
        switch(i) {
            case 0 :
                std::cout << "virtual" << std::endl;
                costmodel = CostModel::VirtualSimulation;
                break;
            case 1 :
                std::cout << "real" << std::endl;
                costmodel = CostModel::RealSimulation;
                break;
            case 2 :
                continue;/*
                std::cout << "real cache flushes" << std::endl;
                costmodel = CostModel::RealSimulationCacheFlushes;*/
                break;
        }
        XFFtCostFactors factors;
        
        std::cout << "no factors" << std::endl;
        analyzeSimulated<T, Allocator, fft::Fastest>(scale*64,
                                                     scale*1984,
                                                     costmodel,
                                                     factors,
                                                     sideEffect);
        
        std::cout << std::endl << "hs 1024 -> 0" << std::endl;
        factors.setMultiplicator(1024,0.f);
        analyzeSimulated<T, Allocator, fft::Fastest>(scale*64,
                                                     scale*1984,
                                                     costmodel,
                                                     factors,
                                                     sideEffect);
    }
    std::cout << "sideEffect " << sideEffect << std::endl;
}
template<typename T>
void printCosts(int const sz) {
    using namespace fft;
    using Tag = Fastest;
    using RealSignalCosts = RealSignalCosts<Tag, T>;
    using RealFBinsCosts = RealFBinsCosts<Tag, T>;
    using AlgoCosts = AlgoCosts<Tag, T>;
    
    int const fft_length = 2*sz;
    
    std::cout << "- for size " << sz << std::endl;

    std::cout << "RealSignalCosts::cost_add_scalar_multiply " << RealSignalCosts::cost_add_scalar_multiply(sz) << std::endl;

    std::cout << "RealSignalCosts::cost_copy " << RealSignalCosts::cost_copy(sz) << std::endl;

    std::cout << "RealSignalCosts::cost_zero_n_raw " << RealSignalCosts::cost_zero_n_raw(sz) << std::endl;

    std::cout << "RealFBinsCosts::cost_mult_assign " << RealFBinsCosts::cost_mult_assign(sz) << std::endl;

    std::cout << "RealFBinsCosts::cost_multiply " << RealFBinsCosts::cost_multiply(sz) << std::endl;

    std::cout << "RealFBinsCosts::cost_multiply_add " << RealFBinsCosts::cost_multiply_add(sz) << std::endl;

    std::cout << "AlgoCosts::cost_fft_inverse " << AlgoCosts::cost_fft_inverse(fft_length) << std::endl;

    std::cout << "AlgoCosts::cost_fft_forward " << AlgoCosts::cost_fft_forward(fft_length) << std::endl;
}
} // NS imajuscule

TEST(BenchmarkConvolutionsScaling, iterateScales_findCheapest) {
    using namespace imajuscule;
    /*
    printCosts<double>(64);
    printCosts<double>(128);
    printCosts<double>(256);
    printCosts<double>(512);
    printCosts<double>(1024);
     */
    smallTest<double, monotonic::aP::Alloc>();
    smallTest<double, a64::Alloc>();
    findCheapest<float, a64::Alloc>();
    findCheapest<double, a64::Alloc>();
}

namespace imajuscule {
std::ostream & operator << (std::ostream & os, const Costs & c) {
    os << "s:" << c.simulated << "\tr:" << c.real << "\tf:" << c.real_flushes;
    return os;
}
}
