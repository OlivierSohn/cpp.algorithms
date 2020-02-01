
namespace imajuscule {

struct Scaling {
    Scaling(int const sz, int const nRepeat) :
    sz(sz), nRepeat(nRepeat) {}
    
    int sz;
    int nRepeat;
    
    bool operator == (Scaling const & o) const {
        return sz==o.sz && nRepeat == o.nRepeat;
    }
    
    bool operator < (Scaling const & o) const {
        if(sz < o.sz) {
            return true;
        }
        if(nRepeat < o.nRepeat) {
            return true;
        }
        return false;
    }

    friend std::ostream & operator << (std::ostream & os, const std::pair<int, Scaling>& p) {
        auto const & s = p.second;
        auto firstSz = p.first;
        auto szMultiplicator = s.sz/firstSz;
        Assert(szMultiplicator * firstSz == s.sz);
        for(int i=0; i<s.nRepeat; ++i) {
            auto const l = szMultiplicator;
            if(l == 1) {
                os << ".";
            }
            else {
                os << "(";
                int nBlanks = l-2;
                if(nBlanks > 0) {
                    os << std::string(nBlanks, ' ');
                }
                os << ")";
            }
        }
        return os;
    }
};

struct ScalingsIterator {
    ScalingsIterator(int firstSz,
                     int totalNCoeffs,
                     std::optional<int> lastSz = {})
    : firstSz(firstSz)
    , totalNCoeffs(totalNCoeffs)
    , lastSz(lastSz)
    {
        if(lastSz) {
            Assert(firstSz <= *lastSz);
            if(firstSz) {
                Assert(firstSz*(*lastSz/firstSz) == *lastSz);
                Assert(is_power_of_two(*lastSz/firstSz));
            }
        }
    }
    
    template<typename F>
    void forEachScaling(F f) {
        if(firstSz <= 0) {
            return;
        }
        std::vector<Scaling> candidate;
        forEachScalingRecursive(candidate, 0, f);
    }

    int const firstSz;
    int const totalNCoeffs;
    std::optional<int> lastSz;

 private:
    template<typename F>
    void forEachScalingRecursive(std::vector<Scaling> pre,
                                 int const preNCoeffs,
                                 F f) {
        if(preNCoeffs >= totalNCoeffs) {
            if(lastSz) {
                if(pre.empty()) {
                    return;
                }
                Assert(pre.back().sz <= *lastSz);
                if(pre.back().sz != *lastSz) {
                    return;
                }
            }
            f(pre);
        }
        else {
            for(auto option : computeNextOptions(pre, preNCoeffs)) {
                auto next = pre;
                auto nextNCoeffs = concatenate(next,
                                               option,
                                               preNCoeffs);
                forEachScalingRecursive(next,
                                        nextNCoeffs,
                                        f);
            }
        }
    }
    
    std::vector<Scaling> computeNextOptions(std::vector<Scaling> const & pre,
                                            int const preNCoeffs)
    {
        Assert(preNCoeffs >= 0);
        Assert(firstSz*(preNCoeffs/firstSz) == preNCoeffs);
        auto const count = (preNCoeffs/firstSz)+1;
        Assert(is_power_of_two(count));
        
        int const remainingCoeffs = totalNCoeffs - preNCoeffs;
        Assert(remainingCoeffs > 0);
        std::vector<Scaling> res;
        if(pre.empty()) {
            Assert(preNCoeffs == 0);
            // add the only valid option here:
            res.emplace_back(firstSz, 1);
        }
        else {
            int const candidateNextSize = count * firstSz;
            // count >= 2 and count is a power of 2, hence candidateNextSize is divisible by 2
            Assert(2*(candidateNextSize/2) == candidateNextSize);

            // If there are enough remaining coefficients, add the option to switch to a bigger size.
            //
            // for example if we have:
            //
            //   1 1 1
            //
            // it means that changing to size "2" at position 2:
            //
            //   1 2
            //
            // was less efficient than repeating twice the size "1".
            //
            // Hence, changing to size "4" at position 4
            //
            //   1 1 1 4
            //
            // will be less efficient than repeating twice the size 1:
            //
            //   1 1 1 1 1
            //
            // but MIGHT be MORE efficient than repeating 3 times the size 1:
            //
            //   1 1 1 1 1 1

            bool const changeSizeIsInteresting = (remainingCoeffs > candidateNextSize/2);

            bool canRepeat = true;
            bool canChangeSize = changeSizeIsInteresting;
            if(lastSz) {
                if(candidateNextSize == *lastSz) {
                    canChangeSize = true;
                    canRepeat = false;
                }
                else if(candidateNextSize > *lastSz) {
                    canChangeSize = false;
                }
            }

            if(canChangeSize) {
                res.emplace_back(candidateNextSize, 1);
            }
            if(canRepeat) {
                Assert(preNCoeffs >= firstSz);
                Assert(count >= 2);

                int const maxNAdditionalRepeats = 1 + ((remainingCoeffs-1) / pre.back().sz);
                Assert(maxNAdditionalRepeats >= 1);
                
                Assert(pre.back().nRepeat >= 1);
                // valid repeats (that allow a change of size afterwards) are 2^n-1: 0, 1, 3, 7, ...
                Assert(is_power_of_two(pre.back().nRepeat+1));

                int const standardNAdditionalRepeats = pre.back().nRepeat + 1;
                Assert(standardNAdditionalRepeats >= 2);
                int const additionalRepeats = std::min(maxNAdditionalRepeats,
                                                       standardNAdditionalRepeats);
                Assert(additionalRepeats >= 1);
                // add the "repeat more of the last size" option:
                res.emplace_back(pre.back().sz, additionalRepeats);
            }
        }
        return res;
    }
    
    static int concatenate(std::vector<Scaling> & v,
                           Scaling const & s,
                           int const preSz)
    {
        Assert(s.nRepeat > 0);
        
        if(!v.empty() && v.back().sz == s.sz) {
            v.back().nRepeat += s.nRepeat;
        }
        else {
            v.push_back(s);
        }
        return preSz + (s.sz * s.nRepeat);
    }
    
};
    
    template<typename ScalingParam>
    auto scalingsToParams(std::vector<Scaling> const & v)
    {
        std::vector<ScalingParam> scalingParams;
        scalingParams.reserve(v.size());
        for(auto const & e:v) {
            ScalingParam s{
                e.sz * e.nRepeat /* count coeffs */,
                {
                    e.sz /* partition size */,
                    e.nRepeat
                }
            };
            scalingParams.push_back(s);
        }
        return scalingParams;
    }

    template<typename C>
    auto mkSimulation(std::vector<Scaling> const & v,
                      int64_t const nCoeffs) {
        using Simulation = typename C::Simulation;
        using ScalingParam = typename Simulation::SetupParam::ScalingParam;
        
        auto scalingParams = scalingsToParams<ScalingParam>(v);
        Simulation sim;
        sim.setup({scalingParams});
        sim.setCoefficientsCount(nCoeffs);
        return sim;
    }
    
    template<typename Simu>
    double virtualCostPerSample(Simu & sim,
                                XFFtCostFactors const & xFftCostFactors) {
        int64_t const end = sim.getBiggestScale();
        if(end) {
            return sim.simuBatch(end,
                                 xFftCostFactors)/end;
        }
        else {
            return 0.;
        }
    }
    
    template<typename C>
    double realCostPerSample(C & c, double & sum) {
        int64_t end = c.getBiggestScale();
        constexpr int minSamples = 40000;
        if(end < minSamples) {
            end *= 1 + (minSamples-1)/end;
        }
        Assert(end >= minSamples);
        
        using namespace imajuscule::profiling;
        auto cost = measure_thread_cpu_one([&sum, end, &c](){
            for(int64_t i= 0; i<end; ++i) {
                sum += c.step({});
            }
        }).count() / 1000000.;
        return cost/end;
    }
    
    template<typename C>
    double realCostPerSampleWithCacheFlushes(C & c,
                                             int64_t const cacheFlushPeriodInSamples,
                                             double & sum) {
        Assert(cacheFlushPeriodInSamples > 0);
        
        using namespace imajuscule::profiling;

        int64_t end = c.getBiggestScale();
        constexpr int minSamples = 40000;
        if(end < minSamples) {
            end *= 1 + (minSamples-1)/end;
        }
        Assert(end >= minSamples);
        
        CpuDuration flush_duration{};
        
        auto test_and_flush_duration = measure_thread_cpu_one([&sum, end, &c, cacheFlushPeriodInSamples, &flush_duration](){
            int64_t remainingToCacheFlush = cacheFlushPeriodInSamples;
            for(int64_t i= 0; i<end; ++i) {
                sum += c.step({});
                --remainingToCacheFlush;
                if(unlikely(!remainingToCacheFlush))
                {
                    std::optional<CpuDuration> local_flush_duration;
                    {
                        Timer t(local_flush_duration);
                        remainingToCacheFlush = cacheFlushPeriodInSamples;
                        polluteCache(sum);
                    }
                    if(unlikely(!local_flush_duration)) {
                        throw std::runtime_error("cannot measure flush duration");
                    }
                    flush_duration += *local_flush_duration;
                }
            }
        });
        auto test_duration = test_and_flush_duration - flush_duration;
        auto cost = test_duration.count() / 1000000.;
        return cost/end;
    }
    
    template<typename Key, typename Value>
    struct VectorIteratorSecondLess {
        using Iterator = typename std::map<Key, Value>::const_iterator;
        bool operator()(Iterator const & it1, Iterator const &it2) {
            return it1->second < it2->second;
        }
    };
            
    template<typename Less, typename Key, typename Value>
    auto orderByValue(std::map<Key, Value> const & results,
                      Less less = VectorIteratorSecondLess<Key, Value>()) {
        std::vector<typename std::remove_reference_t<decltype(results)>::const_iterator> byValue;
        byValue.reserve(results.size());
        for(auto it = results.cbegin(), end = results.cend(); it!=end; ++it) {
            byValue.push_back(it);
        }
        std::sort(byValue.begin(), byValue.end(), less);
        return byValue;
    }

    inline void displayScalingNumerically(std::vector<Scaling> const & r, std::ostream& os) {
    
        int sz = 1;
        for(auto const & s : r) {
            Assert(s.sz >= sz);
            while(s.sz > sz) {
                os << "." << "\t";
                sz *= 2;
            }
            Assert(s.sz == sz);
            os << s.nRepeat << "\t";
            sz *= 2;
        }
    
    }
            
    template<typename Iterator>
    void displayScaling(Iterator it,
                        int const firstSz,
                        int const nCoeffs)
    {
            auto const & [r, cost] = *it;
            std::ostringstream os;
            os << justifyRight(8, std::to_string(firstSz)) << " " << justifyRight(9, std::to_string(nCoeffs)) << " ";
            // cost is in seconds per sample.
            // we put it in microseconds per sample to have a more readable output.
            os << cost * 1000000;
            auto costStr = os.str();
            std::cout << costStr << std::string(std::max(1, 50 - static_cast<int>(costStr.size())), ' ') << "\t";
            
            displayScalingNumerically(r, std::cout);
            
            // graphic display:
            
            /*
             for(auto const & s : r) {
             std::cout << std::make_pair(firstSz, s);
             }*/
            std::cout << std::endl;
    }
            
    template<typename MapIterator>
    auto analyzeScalings(int const firstSz,
                         int const nCoeffs,
                         std::vector<MapIterator> const & byCost,
                         int const countDisplay = -1)
    {
        int i=-1;
        for(auto it:byCost) {
            ++i;
            if(countDisplay > 0 && i>=countDisplay) {
                break;
            }
            displayScaling(it, firstSz, nCoeffs);
        }
    }
            
    template<typename Convolution>
    auto getOptimalScalingScheme_ForTotalCpu_ByVirtualSimulation(int const firstSz,
                                                                 int const nCoeffs,
                                                                 XFFtCostFactors const & xFftCostFactors)
            -> std::optional<std::pair<std::vector<Scaling>, double>>
    {
         std::optional<std::pair<std::vector<Scaling>, double>> best;
            
         ScalingsIterator{
            firstSz,
            nCoeffs
         }.forEachScaling([nCoeffs, &best, &xFftCostFactors](auto const & v){
             auto sim = mkSimulation<Convolution>(v, nCoeffs);
             auto virtualCost = virtualCostPerSample(sim, xFftCostFactors);
             if(!best || virtualCost < best->second) {
                 best = {{v, virtualCost}};
             }
         });
         Assert(best); // assert mis par curiosité pour voir dans quel cas on ne trouve pas de solution
         return best;
     }
}
