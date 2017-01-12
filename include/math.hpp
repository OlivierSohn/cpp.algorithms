
namespace imajuscule {

    constexpr bool zero_or_powOf2(int x) {
        return (x & (x-1)) == 0;
    }
    
    static inline bool zero_or_powOf2(size_t x) {
        while(x) {
            if( 1 == x % 2 ) {
                return false;
            }
            x /= 2;
        }
        return true;
    }
    
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

