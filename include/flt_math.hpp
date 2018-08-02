
namespace imajuscule {
    constexpr float FLOAT_EPSILON = 0.000001f;
    constexpr float FLOAT_ONE_MIN_EPSILON = 1.f - FLOAT_EPSILON;

    template<typename T>
    inline bool same_sign_strict(T a, T b) {
        return (a > static_cast<T>(0) && b > static_cast<T>(0)) || (a < static_cast<T>(0) && b < static_cast<T>(0));
    }

    template<typename T>
    [[nodiscard]] constexpr T clamp_ret(T const v, T const v1, T const v2) {
      Assert(v1 <= v2);
      if(v < v1) {
          return v1;
      }
      if(v > v2) {
          return v2;
      }
      return v;
    }

    template< typename T>
    void clamp_inplace(T & v, T const m, T const M) {
        assert(m <= M);
        if(v > M) {
            v = M;
        }
        else if(v < m) {
            v = m;
        }
    }

} // NS imajuscule
