
namespace imajuscule {

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


}
