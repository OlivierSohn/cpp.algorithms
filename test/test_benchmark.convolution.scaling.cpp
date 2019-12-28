namespace imajuscule {

struct Costs
{
    double simulated, real, real_flushes;
    
    Costs operator * (double v) const {
        return {simulated*v, real*v, real_flushes*v};
    }
    friend std::ostream & operator << (std::ostream & os, const Costs & c) {
        os << "s:" << c.simulated << "\tr:" << c.real << "\tf:" << c.real_flushes;
        return os;
    }
};
                   
template<typename T, typename Tag>
void findCheapest2()
{
    double sideEffect{};
    for(int factor = 100; factor <= 10000; factor *=10)
    {
        for(int firstSz = 1; firstSz < 100000; firstSz*=2)
        {
            int const nCoeffs = factor*firstSz;
            if(nCoeffs > 2000000) {
                continue;
            }
            auto coeffs = mkTestCoeffs<T>(nCoeffs);
           
            std::map<std::vector<Scaling>, Costs> results;
            {
                int count = 0;
                ScalingsIterator{
                    firstSz,
                    static_cast<int>(coeffs.size())
                }.forEachScaling([&coeffs, &count, &results, &sideEffect](auto const & v){
                    auto conv = mkConvolution<T, Tag>(v, coeffs);
                    constexpr int cache_flush_period_samples = 128;
                    results.emplace(v,
                                    Costs{
                        virtualCostPerSample(mkSimulation<T, Tag>(v, coeffs.size())),
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
            analyzeScalings(firstSz, nCoeffs, byCost, nDisplay);
         */
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

} // NS imajuscule

TEST(BenchmarkConvolutionsScaling, iterateScales_findCheapest) {
    using namespace imajuscule;
    findCheapest<float>();
    findCheapest<double>();
}
