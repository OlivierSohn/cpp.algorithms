
namespace imajuscule::audio {

struct Latency {
    constexpr explicit Latency(int latency)
    : latency(latency)
    {}
    
    constexpr int toInteger() const {
        return latency;
    }
    
    Latency operator + (Latency const & o) const {
        return Latency(latency + o.latency);
    }
    Latency operator - (Latency const & o) const {
        return Latency(latency - o.latency);
    }
    
    constexpr bool operator == (Latency const & o) const {
        return latency == o.latency;
    }
    constexpr bool operator != (Latency const & o) const {
        return latency != o.latency;
    }
    constexpr bool operator < (Latency const & o) const {
        return latency < o.latency;
    }
    constexpr bool operator > (Latency const & o) const {
        return latency > o.latency;
    }
    constexpr bool operator <= (Latency const & o) const {
        return latency <= o.latency;
    }
    constexpr bool operator >= (Latency const & o) const {
        return latency >= o.latency;
    }
private:
    int latency;
};

using FftSpecs = std::map<uint32_t, int>; // sz -> historySize

struct MinSizeRequirement {
    MinSizeRequirement(int minXSize,
                       int minYSize,
                       FftSpecs xFftSizes,
                       int minWorkSize)
    : minXSize(minXSize)
    , minYSize(minYSize)
    , xFftSizes(xFftSizes)
    , minWorkSize(minWorkSize)
    {}
    
    /*
     The needed size of the x buffer is deduced both from the max size of the ffts, and from this parameter.
     */
    int minXSize;
        
    int minYSize;
    
    FftSpecs xFftSizes;
    
    int minWorkSize;
    
    void mergeWith(MinSizeRequirement const & o) {
        minWorkSize = std::max(minWorkSize,
                               o.minWorkSize);
        minXSize = std::max(minXSize,
                            o.minXSize);
        minYSize = std::max(minYSize,
                            o.minYSize);
        
        for(auto const & [sz, historySize] : o.xFftSizes) {
            auto [it, emplaced] = xFftSizes.try_emplace(sz, historySize);
            if(!emplaced) {
                it->second = std::max(it->second, historySize);
            }
        }
    }
};

template<typename A, typename T, typename FFTTag>
struct Simulation_;

template<typename A, typename T, typename FFTTag>
using Simulation = typename Simulation_<A, T, FFTTag>::type;

}
