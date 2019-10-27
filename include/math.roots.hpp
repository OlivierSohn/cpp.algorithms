/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule {

// Assumes that f(curFrom) and f(curTo) have opposite signs, and that f is continuous.
template<typename F, typename FDER, typename T = decltype(std::declval<F>()({}))>
auto find_root(F f, FDER fder, T curFrom, T curTo) -> std::optional<T> {
    static_assert(std::is_same_v<
                  decltype(std::declval<F>()({})),
                  decltype(std::declval<FDER>()({}))
                  >);
    
    struct XY {
        XY(F f, T x) :
        x(x), y(f(x))
        {}
        
        XY(T x, T y) :
        x(x), y(y)
        {}
        
        T x,y;
    };

    struct Bounds {
        Bounds(F f, T x1, T x2)
        : positive(f, x1)
        , negative(f, x2)
        {}
        
        T middle() const {
            return 0.5 * (positive.x + negative.x);
        }
        T falseposition_mid() const {
            auto denom = positive.y - negative.y;
            assert(denom > 0);
            return (negative.x*positive.y - positive.x*negative.y) / denom;
        }
        
        bool xStrictlyWithin(T x) const {
            if(negative.x < positive.x) {
                return negative.x < x && x < positive.x;
            }
            else {
                return positive.x < x && x < negative.x;
            }
        }

        bool update(XY const & v) {
            if(negative.x < positive.x) {
                if(v.y > 0) {
                    if( v.x < positive.x) {
                        positive = v;
                    }
                    else {
                        return false;
                    }
                }
                else {
                    if( v.x > negative.x) {
                        negative = v;
                    }
                    else {
                        return false;
                    }
                }
            }
            else {
                if(v.y > 0) {
                    if( v.x > positive.x) {
                        positive = v;
                    }
                    else {
                        return false;
                    }
                }
                else {
                    if( v.x < negative.x) {
                        negative = v;
                    }
                    else {
                        return false;
                    }
                }
            }
            return true;
        };
        XY positive, negative;
    } bounds(f, curFrom, curTo);

    auto constexpr epsilon = std::numeric_limits<T>::epsilon();

    if(std::abs(bounds.negative.y) < epsilon) {
        return bounds.negative.x;
    }
    if(std::abs(bounds.positive.y) < epsilon) {
        return bounds.positive.x;
    }
    
    if(bounds.negative.y > 0) {
        if(bounds.positive.y > 0) {
            // invalid bounds
            return {};
        }
        std::swap(bounds.negative, bounds.positive);
    }
    else if(bounds.positive.y < 0) {
        // invalid bounds
        return {};
    }
    
    enum class Method {
        Newton,
        FalsePosition
    } method = Method::FalsePosition;

    auto performNewton = [] (T const x, T const y, T const slope) -> std::optional<T> {
        if(slope) {
            return {x-y/slope};
        }
        return {};
    };
    

    T X, Y;
    while(true) {
        if(likely(method == Method::Newton)) {
            T const newtonSlope = fder(X);
            auto mayNewtonX = performNewton(X, Y, newtonSlope);
            if(!mayNewtonX) {
                // zero slope
                method = Method::FalsePosition;
                continue;
                
            }
            X = *mayNewtonX;
            Y = f(X);
            if(std::abs(Y) < epsilon) {
                return { X };
            }
            if(!bounds.update({X, Y})) {
                // outside bounds
                method = Method::FalsePosition;
                continue;
            }
        }
        else {
            assert(method==Method::FalsePosition);
            X = bounds.falseposition_mid();
            if(!bounds.xStrictlyWithin(X)) {
                return { bounds.middle() };
            }
            Y = f(X);
            if(std::abs(Y) < epsilon) {
                return { X };
            }
            if(!bounds.update({X, Y})) {
                // by design this should never happen
                assert(0);
                return {};
            }
            method = Method::Newton;
        }
    }
    return {};
}

} // NS imajuscule
