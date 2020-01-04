/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule {
  /*
  Considers that the range of unsigned values are on a circle, and measures the min
  distance between values on that circle, hence:

  {
    unsigned int a = 3;
    unsigned int b = 5;
    Assert(cyclic_unsigned_dist(a,b) == cyclic_unsigned_dist(b,a) == 2);
  }
  {
    unsigned int a = std::numeric_limits<unsigned int>::min();
    unsigned int b = std::numeric_limits<unsigned int>::max();
    Assert(cyclic_unsigned_dist(a,b) == cyclic_unsigned_dist(b,a) == 1);
  }
  Assert(1 == cyclic_unsigned_dist(std::numeric_limits<unsigned int>::min(), std::numeric_limits<unsigned int>::max()))
  */
  template<typename U>
  constexpr U cyclic_unsigned_dist(U a, U b) {
    static_assert(std::is_unsigned_v<U>);
    return std::min(a-b, b-a); // it is intended that one of them can underflow.
  }

  constexpr bool is_even(unsigned int i) {
    return (i & 1) == 0;
  }
  constexpr bool is_odd(unsigned int i) {
    return !is_even(i);
  }

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

    static constexpr inline size_t pow2(size_t power) {
        size_t res = (static_cast<size_t>(1) << power);
        return res;
    }

    template<uint64_t power>
    static constexpr inline uint64_t pow2() {
        constexpr uint64_t res = (((uint64_t)1) << power);
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
      return __builtin_ctz(x);
//        return 32 - count_leading_zeroes(~x & (x-1));
    }

static inline int count_set_bits(uint32_t x) {
    return __builtin_popcount(x);
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

template<typename T>
int countDecimalNumbersBeforeTheDot(T val){
    val = std::abs(val);
    int count = 0;
    while(true) {
        if(val < 1) {
            return count;
        }
        val /= 10;
        ++count;
    }
}

static inline int gcd(int a, int b)
{
    if (a == 0) {
        return b;
    }
    return gcd(b % a, a);
}
  
static inline int64_t ppcm(int ia, int ib)
{
    int denom = gcd(ia,ib);
    if(0 == denom) {
        return 0;
    }
    int64_t a = ia;
    int64_t b = ib;
    int64_t prod = a*b;
    return std::abs(prod) / denom;
}

} // NS imajuscule
