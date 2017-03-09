
namespace imajuscule {
    constexpr float FLOAT_EPSILON = 0.000001f;
    constexpr float FLOAT_ONE_MIN_EPSILON = 1.f - FLOAT_EPSILON;

    inline bool same_sign_strict(float a, float b) {
        return (a > 0.f && b > 0.f) || (a < 0.f && b < 0.f);
    }
    
    constexpr float clamp(float v, float v1, float v2) {
        if(v < v1) {
            return v1;
        }
        if(v> v2) {
            return v2;
        }
        return v;
    }
    
} // NS imajuscule

