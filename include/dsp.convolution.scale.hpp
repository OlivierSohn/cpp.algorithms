namespace imajuscule::audio {

struct CountDroppedScales {
    constexpr explicit CountDroppedScales(int n)
    : n(n)
    {}
    
    constexpr int toInteger() const {
        return n;
    }
private:
    int n;
};

struct ScaleConvolution_ {
    static constexpr int latencyForDroppedConvolutions(CountDroppedScales const & nDropped) {
        return static_cast<int>(pow2(static_cast<size_t>(nDropped.toInteger())))-1;
    }
    
    /*
     * The default value is optimal for the system I developped on (OSX / Macbook Air 2015 / intel core i7).
     * see 'TEST(BenchmarkConvolutions, scaled_vs_brute)' on how to determine the best value for a given system.
     *
     * TODO ideally this value should be a global, computed at initialization time.
     */
    static constexpr CountDroppedScales nDroppedOptimalFor_Split_Bruteforce_Fft = CountDroppedScales(6);
};


constexpr int countPartitions(int nCoeffs, int partition_size) {
    if(nCoeffs <= 0) {
        return 0;
    }
    Assert(partition_size);
    return 1 + (nCoeffs-1)/partition_size;
}

namespace SameSizeScales {
    /*
     Assuming that scales have resolutions of 1,2,4,8...
     for every scale, the number of blocks = { N(k), k in 1 .. S }
     Assuming that between scale i and (i+1) there are nOverlap*2^(i-1) blocks in common
     N = sum (n in 0 .. (S-1)) 2^n N(n+1) - scaleFadeSz::inSmallerUnits * sum (n in 1.. (S-1)) 2^(n-1)
     and since we want the number of blocks to be equal
     (hence the namespace name 'SameSizeScales'), we will solve this:
     N = NF * sum (n in 0 .. (S-1)) 2^n - scaleFadeSz::inSmallerUnits * sum (n in 0.. (S-2)) 2^n
     N = NF * (2^S - 1) - scaleFadeSz::inSmallerUnits * (2^(S-1)-1)
     NF = (N + scaleFadeSz::inSmallerUnits * (2^(S-1)-1)) / (2^S - 1)
     */
static inline int get_scale_sz(int response_sz, int n_scales) {
    assert(response_sz>=0);
    if(n_scales <= 1) {
        return response_sz;
    }
    int numerator = response_sz + static_cast<int>(scaleFadeSz::inSmallerUnits * (static_cast<int>(pow2(n_scales-1)) - 1));
    int denominator = static_cast<int>(pow2(n_scales)) - 1;
    int res = countPartitions(numerator, denominator);
    if(res%2) {
        return res + 1;
    }
    return res;
}
static inline int get_max_response_sz(int n_scales, int sz_one_scale) {
    return
    sz_one_scale * (pow2(n_scales) - 1) -
    scaleFadeSz::inSmallerUnits * (static_cast<int>(pow2(n_scales-1)) - 1);
}

int constexpr getDelays(int scale_sz, Latency latency) {
    // split = nFadeIn + latB - latA
    // scale_sz = nFadeIn + (1 + 2*(delay + latA)) - latA
    // scale_sz - 1 - nFadeIn = 2*delay + latA
    // delay = 0.5 * (scale_sz - 1 - nFadeIn - latA)
    assert(0 == (scale_sz-scaleFadeSz::inSmallerUnits) % 2);
    return (scale_sz - 1 - scaleFadeSz::inSmallerUnits - latency.toInteger()) / 2;
}
};

}
