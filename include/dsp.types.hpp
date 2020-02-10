
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

/*
 For an intuitive presentation of the 2 methods, see:
 https://www.slideshare.net/GourabGhosh4/overlap-add-overlap-savedigital-signal-processing
 
 For our use-case, overlap save has several advantages:
 . it allows to perform ffts of x _concurrently_ using a single x buffer because we don't need 0-padding
 . it makes less (half) memory writes to y (because only the second half of the spectrum is used)
 . it uses less (half) memory for y (for the same reason mentionned above)
 
 The only inconvenient overlap save has is that:
 . with overlap add, we can do the x 0-padding incrementally, so the worst (peak) number of reads / writes to x is
   maxFftSize/4
 . with overlap save, copying x needs to happen at once, so the worst (peak) number of reads / writes to x is
   maxFftSize/2
 
 When the biggest fft of x and the write of the results to y happen in the same step,
 this inconvenient is entirely compensated by the fact that
 there are less memory writes to y.
 
 But when the biggest fft of x and the write to y happen in _different_ steps
 (like in the finegrained case), this inconvenient leads to a higher worst case step time.
*/
enum class Overlap {
    Add,
    Save
};


}
