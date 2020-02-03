
namespace imajuscule {


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
    MinSizeRequirement() = delete;
    
    /*
     The needed size of the x buffer is deduced both from the max size of the ffts, and from this parameter.
     */
    int minXSize;
        
    int minYSize;
    // "Anticipated" writes touch the block (of size minYSize) after the current block.
    int minYAnticipatedWrites;
    
    FftSpecs xFftSizes;
    
    void mergeWith(MinSizeRequirement const & o) {
        minXSize = std::max(minXSize, o.minXSize);
        
        if(!minYSize || !o.minYSize) {
            minYSize = std::max(minYSize, o.minYSize);
        }
        else
        {
            int res = static_cast<int>(ppcm(minYSize,
                                            o.minYSize));
            
            // only because in practice we have powers of 2.
            Assert(std::max(minYSize,
                            o.minYSize) == res);

            minYSize = res;
        }
        
        minYAnticipatedWrites = std::max(minYAnticipatedWrites,
                                         o.minYAnticipatedWrites);
        
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
