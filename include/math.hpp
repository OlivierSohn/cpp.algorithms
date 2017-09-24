/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

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
    constexpr void increment(E & e) {
        e = static_cast<E>(to_underlying(e) + 1);
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
    
    template<int64_t power>
    static constexpr int64_t pow2() {
        constexpr int64_t res = (((int64_t)1) << power);
        return res;
    }
    
    static inline uint32_t count_leading_zeroes(uint32_t arg) {
        if (arg == 0) {
            return 32;
        }
        return __builtin_clz(arg);
    }
    
    static inline uint32_t floor_power_of_two(uint32_t arg) {
        return arg ? ceil_power_of_two(1+arg/2) : 0;
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
    
    template<typename T>
    T exp_mean(T a, T b) {
        if(a < 0 && b < 0) {
            return -exp_mean(-a, -b);
        }
        if(a < 0 || b < 0) {
            throw std::logic_error("logarithmic mean with negative input");
        }
        if(a == b) {
            return a;
        }
        float mid_log = (std::log(static_cast<float>(b)) + std::log(static_cast<float>(a))) / 2.f;
        T res = std::exp(mid_log);
        
        return static_cast<T>(res);
    }
    
    template<typename ITER, typename VAL = typename ITER::value_type>
    VAL dichotomic_sum(ITER it, ITER end) {
        auto d = std::distance(it, end);
        if(0 == d) {
            return {};
        }
        if(1 == d) {
            return *it;
        }
        auto mid = it + d/2;
        return dichotomic_sum(it, mid) + dichotomic_sum(mid, end);
    }
    
    template< typename T>
    void clamp(T & v, T m, T M) {
        assert(m <= M);
        if(v > M) {
            v = M;
        }
        else if(v < m) {
            v = m;
        }
    }

} // NS imajuscule

