/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule
{
    struct itp
    {
        // WARNING if you modify this enum, modify also the floats below
        enum interpolation : unsigned char
        {
            INTERPOLATION_LOWER_BOUND = 0,

            LINEAR = INTERPOLATION_LOWER_BOUND,
            PROPORTIONAL_VALUE_DERIVATIVE,
            EASE_IN_first_value,
            EASE_IN_QUAD,
            EASE_IN_CUBIC,
            EASE_IN_QUART,
            EASE_IN_QUINT,
            EASE_IN_SINE,
            EASE_IN_EXPO,
            EASE_IN_CIRC,
            EASE_IN_last_value,
            EASE_OUT_first_value,
            EASE_OUT_QUAD,
            EASE_OUT_CUBIC,
            EASE_OUT_QUART,
            EASE_OUT_QUINT,
            EASE_OUT_SINE,
            EASE_OUT_EXPO,
            EASE_OUT_CIRC,
            EASE_OUT_last_value,
            EASE_INOUT_first_value,
            EASE_INOUT_QUAD,
            EASE_INOUT_CUBIC,
            EASE_INOUT_QUART,
            EASE_INOUT_QUINT,
            EASE_INOUT_SINE,
            EASE_INOUT_EXPO,
            EASE_INOUT_CIRC,
            EASE_INOUT_last_value,

            INTERPOLATION_UPPER_BOUND,
        };
        
        static const float NO_EASE;
        static const float EASE_IN;
        static const float EASE_OUT;
        static const float EASE_IO;
        
        // Naming convention: polynomial is "P{degree}"

        // values for NO_EASE:
        static const float P1;
        static const float PVD;

        // values for EASE_IN/OUT/IO:
        static const float P2;
        static const float P3;
        static const float P4;
        static const float P5;
        static const float SINE;
        static const float EXPO;
        static const float CIRC;
        
        
        static constexpr bool isValid( interpolation type )
        {
            switch(type)
            {
                default:
                    return false;
                case LINEAR:
                    return true;
                case PROPORTIONAL_VALUE_DERIVATIVE:
                    return true;
                case EASE_IN_QUAD:
                    return true;
                case EASE_IN_CUBIC:
                    return true;
                case EASE_IN_QUART:
                    return true;
                case EASE_IN_QUINT:
                    return true;
                case EASE_IN_SINE:
                    return true;
                case EASE_IN_EXPO:
                    return true;
                case EASE_IN_CIRC:
                    return true;
                case EASE_OUT_QUAD:
                    return true;
                case EASE_OUT_CUBIC:
                    return true;
                case EASE_OUT_QUART:
                    return true;
                case EASE_OUT_QUINT:
                    return true;
                case EASE_OUT_SINE:
                    return true;
                case EASE_OUT_EXPO:
                    return true;
                case EASE_OUT_CIRC:
                    return true;
                case EASE_INOUT_QUAD:
                    return true;
                case EASE_INOUT_CUBIC:
                    return true;
                case EASE_INOUT_QUART:
                    return true;
                case EASE_INOUT_QUINT:
                    return true;
                case EASE_INOUT_SINE:
                    return true;
                case EASE_INOUT_EXPO:
                    return true;
                case EASE_INOUT_CIRC:
                    return true;
            }
            return false;
        }
        
        static const enumTraversal & interpolation_traversal();

        // returns true if value is an interpolation
        static bool intIsReal(int val);
        static interpolation intToInterpolation(int val, bool & bReal);
        static const char * interpolationInfo(int);

        // cf. http://gizma.com/easing/
        
        // parameters are :
        // t : time (from 0 to d)
        // d : duration
        // b : value at t==0
        // c : added value i.e b+c is the value at t==d
        
        // wraps the other functions
        static float interpolate(interpolation type, float t_start, float t_cur, float t_end, float val_start, float val_end) {
            return interpolate(type
                               , t_cur - t_start
                               , val_start
                               , val_end - val_start
                               , t_end - t_start);
        }
        
        static float interpolate(interpolation type, float t_start_to_cur, float val_start = 0.f, float val_range = 1.f, float t_start_to_end = 1.f);
       
        static float linearTween(float t, float b, float c, float d);
        static float proportional_value_derivative(float t, float b, float c, float d);

        // quadratic easing in - accelerating from zero velocity
        
        static float easeInQuad(float t, float b, float c, float d);
        
        // quadratic easing out - decelerating to zero velocity
        
        static float easeOutQuad(float t, float b, float c, float d);
        
        // quadratic easing in/out - acceleration until halfway, then deceleration
        
        static float easeInOutQuad(float t, float b, float c, float d);
        
        // cubic easing in - accelerating from zero velocity
        
        static float easeInCubic(float t, float b, float c, float d);
        
        // cubic easing out - decelerating to zero velocity
        
        static float easeOutCubic(float t, float b, float c, float d);
        
        // cubic easing in/out - acceleration until halfway, then deceleration
        
        static float easeInOutCubic(float t, float b, float c, float d);
        
        // quartic easing in - accelerating from zero velocity
        
        static float easeInQuart(float t, float b, float c, float d);
        
        // quartic easing out - decelerating to zero velocity
        
        static float easeOutQuart(float t, float b, float c, float d);
        
        // quartic easing in/out - acceleration until halfway, then deceleration
        
        static float easeInOutQuart(float t, float b, float c, float d);
        
        // quintic easing in - accelerating from zero velocity
        
        static float easeInQuint(float t, float b, float c, float d);
        
        // quintic easing out - decelerating to zero velocity
        
        static float easeOutQuint(float t, float b, float c, float d);
        
        // quintic easing in/out - acceleration until halfway, then deceleration
        
        static float easeInOutQuint(float t, float b, float c, float d);
        
        // sinusoidal easing in - accelerating from zero velocity
        
        static float easeInSine(float t, float b, float c, float d);
        
        // sinusoidal easing out - decelerating to zero velocity
        
        static float easeOutSine(float t, float b, float c, float d);
        
        // sinusoidal easing in/out - accelerating until halfway, then decelerating
        
        static float easeInOutSine(float t, float b, float c, float d);
        
        // exponential easing in - accelerating from zero velocity
        
        static float easeInExpo(float t, float b, float c, float d);
        
        // exponential easing out - decelerating to zero velocity
        
        static float easeOutExpo(float t, float b, float c, float d);
        
        // exponential easing in/out - accelerating until halfway, then decelerating
        
        static float easeInOutExpo(float t, float b, float c, float d);
        
        // circular easing in - accelerating from zero velocity
        
        static float easeInCirc(float t, float b, float c, float d);
        
        // circular easing out - decelerating to zero velocit
        
        static float easeOutCirc(float t, float b, float c, float d);
        
        // circular easing in/out - acceleration until halfway, then deceleration
        
        static float easeInOutCirc(float t, float b, float c, float d);

        private:
            static enumTraversal gInterpolationTraversal;
    };
    
    template<typename T = float>
    struct NormalizedInterpolation {
        using Tr = NumTraits<T>;
        static constexpr T zero = Tr::zero();
        
        NormalizedInterpolation() = default;
        
        NormalizedInterpolation(itp::interpolation i) {
            setInterpolation(i);
        }
        
        void setInterpolation(itp::interpolation i) {
            interp = i;
            assert(itp::intIsReal(interp));
        }
        
        T get_value(T t_cur, T t_end, T val_start, T val_end) const {
            assert(t_cur >= zero); // else need to return val_start
            assert(itp::intIsReal(interp));
            
            if(t_cur > t_end) {
                return val_end;
            }
            
            return get_unfiltered_value(t_cur,
                                        t_end,
                                        val_start,
                                        val_end);
        }
        
        T get_unfiltered_value(T t_cur, T t_end, T val_start, T val_end) const {
            return itp::interpolate(interp,
                                    zero /* t start */,
                                    t_cur,
                                    t_end,
                                    val_start,
                                    val_end);
        }
        
    private:
        itp::interpolation interp : 5;
    };
}
