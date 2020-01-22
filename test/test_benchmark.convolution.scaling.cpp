namespace imajuscule {

struct Costs
{
    double simulated, real, real_flushes;
    
    Costs operator * (double v) const {
        return {simulated*v, real*v, real_flushes*v};
    }
};

template<typename T, typename Tag>
void findCheapest2(int const firstSz, int const nCoeffs, XFFtCostFactors const & factors, double & sideEffect)
{
    using C = CustomScaleConvolution<FFTConvolutionIntermediate < PartitionnedFFTConvolutionCRTP<T, Tag> >>;
    
    auto coeffs = mkTestCoeffs<T>(nCoeffs);
    
    std::map<std::vector<Scaling>, Costs> results;
    {
        int count = 0;
        ScalingsIterator{
            firstSz,
            static_cast<int>(coeffs.size())
        }.forEachScaling([&coeffs, &count, &results, &factors, &sideEffect](auto const & v){
            auto conv = mkConvolution<T, Tag>(v, coeffs);
            constexpr int cache_flush_period_samples = 128;
            auto sim = mkSimulation<C>(v, coeffs.size());
            results.emplace(v,
                            Costs{
                virtualCostPerSample(sim, factors),
                realCostPerSample(conv, sideEffect),
                realCostPerSampleWithCacheFlushes(conv, cache_flush_period_samples, sideEffect)
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
template<typename T, typename Tag>
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
            findCheapest2<T, Tag>(firstSz, nCoeffs, factors, sideEffect);
        }
    }
    std::cout << "sideEffect " << sideEffect << std::endl;
}


template<typename T>
void findCheapest()
{
    for_each(fft::Tags, [](auto t) {
        using Tag = decltype(t);
        COUT_TYPE(Tag); std::cout << std::endl;
        findCheapest2<T, Tag>();
    });
}
enum class CostModel {
    VirtualSimulation,
    RealSimulation,
    RealSimulationCacheFlushes
};

template<typename T, typename Tag>
void analyzeSimulated(int const firstSz, int const nCoeffs, CostModel model, XFFtCostFactors const & factors, double & sideEffect)
{
    using C = CustomScaleConvolution<FFTConvolutionIntermediate < PartitionnedFFTConvolutionCRTP<T, Tag> >>;
    
    auto coeffs = mkTestCoeffs<T>(nCoeffs);
    
    std::map<std::vector<Scaling>, double> results;
    {
        int count = 0;
        ScalingsIterator{
            firstSz,
            static_cast<int>(coeffs.size())
        }.forEachScaling([model, &coeffs, &count, &results, &factors, &sideEffect](auto const & v){
            displayScalingNumerically(v, std::cout);
            std::cout << std::endl;
            if(v.size() != 1) {
                return;
            }
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
                    C conv = mkConvolution<T, Tag>(v, coeffs);
                    results.emplace(v,
                                    realCostPerSample(conv, sideEffect));
                    
                }
                    break;
                case CostModel::RealSimulationCacheFlushes:
                {
                    C conv = mkConvolution<T, Tag>(v, coeffs);
                    constexpr int cache_flush_period_samples = 128;
                    results.emplace(v,
                                    realCostPerSampleWithCacheFlushes(conv, cache_flush_period_samples, sideEffect));
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
                std::cout << "real cache flushes" << std::endl;
                costmodel = CostModel::RealSimulationCacheFlushes;
                break;
        }
        XFFtCostFactors factors;
        
        std::cout << "no factors" << std::endl;
        analyzeSimulated<double, fft::Fastest>(scale*64,
                                               scale*1984,
                                               costmodel,
                                               factors,
                                               sideEffect);
        
        /*            std::cout << std::endl << "hs 1024 -> 0" << std::endl;
         factors.setMultiplicator(1024,0.f);
         analyzeSimulated<double, fft::Fastest>(64, 1984, costmodel, factors, sideEffect);
         */
    }
    std::cout << "sideEffect " << sideEffect << std::endl;
}
} // NS imajuscule

TEST(BenchmarkConvolutionsScaling, iterateScales_findCheapest) {
    using namespace imajuscule;
    smallTest();
    findCheapest<float>();
    findCheapest<double>();
}

namespace imajuscule {
std::ostream & operator << (std::ostream & os, const Costs & c) {
    os << "s:" << c.simulated << "\tr:" << c.real << "\tf:" << c.real_flushes;
    return os;
}
}
