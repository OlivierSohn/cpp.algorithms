
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
    
    constexpr unsigned int ceil_power_of_two(unsigned int v)
    {
        --v;
        v |= v >> 1;
        v |= v >> 2;
        v |= v >> 4;
        v |= v >> 8;
        v |= v >> 16;
        v++;
        return v;
    }

    template <typename E>
    constexpr typename std::underlying_type<E>::type to_underlying(E e) {
        return static_cast<typename std::underlying_type<E>::type>(e);
    }
    
    template <typename E>
    constexpr int min_bits() {
        return ceil_power_of_two(to_underlying(E::SIZE_ENUM));
    }
    
    template<int N>
    constexpr uintptr_t removeLowBits(uintptr_t p) {
        return (p >> N) << N;
    }
    
    constexpr bool is_power_of_two(size_t n) { return ((n != 0) && !(n & (n - 1))); }
    
    static constexpr size_t pow2(int power) {
        size_t res = (1 << power);
        return res;
    }
    
    static inline uint32_t count_leading_zeroes(uint32_t arg) {
        if (arg == 0) {
            return 32;
        }
        return __builtin_clz(arg);
    }
    
    // count trailing zeroes
    static inline uint32_t count_trailing_zeroes(uint32_t x)
    {
        return 32 - count_leading_zeroes(~x & (x-1));
    }
    
    template <typename T>
    constexpr T expt_unsigned(T p, unsigned int q) {
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
    
    template <int q, typename T>
    T expt(T p) {
        static_assert(!std::is_integral<T>(), "use expt_unsigned instead");
        if(q == 1) {
            return p;
        }
        if(q == 0) {
            return NumTraits<T>::one();
        }
        if(q >= 0) {
            return expt_unsigned(p, static_cast<unsigned int>(q));
        }
        return expt_unsigned(1/p, static_cast<unsigned int>(-q));
    }
    
} // NS imajuscule

