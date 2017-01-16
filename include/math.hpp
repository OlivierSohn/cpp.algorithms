
namespace imajuscule {

    constexpr unsigned int relevantBits(unsigned int v) {
        unsigned int n = 0;
        while(v) {
            ++n;
            v >>= 1;
        }
        return n;
    }
    
    constexpr unsigned int power_of_two_exponent(unsigned int v) {
        unsigned int n = 0;
        while (v >>= 1) {
            ++n;
        }
        return n;
    }
    
    constexpr bool zero_or_powOf2(int x) {
        return (x & (x-1)) == 0;
    }
    
    constexpr bool is_power_of_two(size_t n) { return ((n != 0) && !(n & (n - 1))); }
    
    static constexpr size_t pow2(int power) {
        size_t res = (1 << power);
        return res;
    }
    
    template <typename T>
    T expt_unsigned(T p, unsigned int q) {
        assert(q < 100000); // else probably a sign error.
        T r(1);
        
        while (q != 0) {
            if (q % 2 == 1) {
                r *= p;
                q--;
            }
            p *= p;
            q /= 2;
        }
        
        return r;
    }
    
    template <typename T>
    T expt(T p, int q) {
        static_assert(!std::is_integral<T>(), "use expt_unsigned instead");
        if(q >= 0) {
            return expt_unsigned(p, static_cast<unsigned int>(q));
        }
        return expt_unsigned(1/p, static_cast<unsigned int>(-q));
    }
    
} // NS imajuscule

