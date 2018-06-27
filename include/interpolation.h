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
      EASE_IN_QUAD, // = 3
      EASE_IN_CUBIC,
      EASE_IN_QUART,
      EASE_IN_QUINT,
      EASE_IN_SINE,
      EASE_IN_EXPO,
      EASE_IN_CIRC,
      EASE_IN_last_value,
      EASE_OUT_first_value,
      EASE_OUT_QUAD, // = 12
      EASE_OUT_CUBIC,
      EASE_OUT_QUART,
      EASE_OUT_QUINT,
      EASE_OUT_SINE,
      EASE_OUT_EXPO,
      EASE_OUT_CIRC,
      EASE_OUT_last_value,
      EASE_INOUT_first_value,
      EASE_INOUT_QUAD, // = 21
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


    static constexpr bool intIsReal(int val)
    {
      switch (static_cast<itp::interpolation> (val))
      {
      case EASE_IN_first_value:
      case EASE_IN_last_value:
      case EASE_OUT_first_value:
      case EASE_OUT_last_value:
      case EASE_INOUT_first_value:
      case EASE_INOUT_last_value:
        return false;
      default:
        return (val >= INTERPOLATION_LOWER_BOUND) && (val < INTERPOLATION_UPPER_BOUND);
      }
    }

    // returns true if value is an interpolation
    static interpolation toItp(int i) {
      auto v = static_cast<interpolation>(i);
      if (intIsReal(v)) {
        return v;
      }
      Assert(0 && "interpolation is not valid");
      return itp::LINEAR;
    }

    static const enumTraversal & interpolation_traversal();

    static interpolation intToInterpolation(int val, bool & bReal);
    static const char * interpolationInfo(int);

    // cf. http://gizma.com/easing/

    // parameters are :
    // t : time (from 0 to d)
    // d : duration
    // b : value at t==0
    // c : added value i.e b+c is the value at t==d

    // wraps the other functions
    template<typename T>
    static T interpolate(interpolation type, T t_start, T t_cur, T t_end, T val_start, T val_end) {
      return interpolate(type
                         , t_cur - t_start
                         , val_start
                         , val_end - val_start
                         , t_end - t_start);
    }

    template<typename T>
    static T interpolate( interpolation type, T t, T b, T c, T d)
    {
      switch(type)
      {
        case LINEAR:
          return linearTween(t,b,c,d);
        case PROPORTIONAL_VALUE_DERIVATIVE:
          return proportional_value_derivative(t,b,c,d);
        case EASE_IN_QUAD:
          return easeInOrd2(t,b,c,d);
        case EASE_IN_CUBIC:
          return easeInOrd3(t,b,c,d);
        case EASE_IN_QUART:
          return easeInOrd4(t,b,c,d);
        case EASE_IN_QUINT:
          return easeInOrd5(t,b,c,d);
        case EASE_IN_SINE:
          return easeInSine(t,b,c,d);
        case EASE_IN_EXPO:
          return easeInExpo(t,b,c,d);
        case EASE_IN_CIRC:
          return easeInCirc(t,b,c,d);
        case EASE_OUT_QUAD:
          return easeOutOrd2(t,b,c,d);
        case EASE_OUT_CUBIC:
          return easeOutOrd3(t,b,c,d);
        case EASE_OUT_QUART:
          return easeOutOrd4(t,b,c,d);
        case EASE_OUT_QUINT:
          return easeOutOrd5(t,b,c,d);
        case EASE_OUT_SINE:
          return easeOutSine(t,b,c,d);
        case EASE_OUT_EXPO:
          return easeOutExpo(t,b,c,d);
        case EASE_OUT_CIRC:
          return easeOutCirc(t,b,c,d);
        case EASE_INOUT_QUAD:
          return easeInOutOrd2(t,b,c,d);
        case EASE_INOUT_CUBIC:
          return easeInOutOrd3(t,b,c,d);
        case EASE_INOUT_QUART:
          return easeInOutOrd4(t,b,c,d);
        case EASE_INOUT_QUINT:
          return easeInOutOrd5(t,b,c,d);
        case EASE_INOUT_SINE:
          return easeInOutSine(t,b,c,d);
        case EASE_INOUT_EXPO:
          return easeInOutExpo(t,b,c,d);
        case EASE_INOUT_CIRC:
          return easeInOutCirc(t,b,c,d);
        default:
          assert(0);
          return static_cast<T>(0);
      }
    }
    //simple linear tweening - no easing, no acceleration


    template<typename T>
    static T linearTween(
        T t, // t_cur - t_start
        T b, // val_start
        T c, // val_end - val_start
        T d  // t_end - t_start
        )
    {
      return c*t/d + b;
    }


    // quadratic easing In - accelerating from zero velocity


    template<typename T>
    static T easeInOrd2(T t, T b, T c, T d) {
      t /= d;
      return c*t*t + b;
    }


    // quadratic easing Out - decelerating to zero velocity

    template<typename T>
    static T easeOutOrd2(T t, T b, T c, T d) {
      t /= d;
      return -c * t*(t-static_cast<T>(2)) + b;
    }



    // quadratic easing In/Out - acceleration until halfway, then deceleration

    template<typename T>
    static T easeInOutOrd2(T t, T b, T c, T d) {
      t /= d/static_cast<T>(2);
      if (t < static_cast<T>(1)) return c/static_cast<T>(2)*t*t + b;
      t--;
      return -c/static_cast<T>(2) * (t*(t-static_cast<T>(2)) - static_cast<T>(1)) + b;
    }


    // Ord3 easing In - accelerating from zero velocity

    template<typename T>
    static T easeInOrd3(T t, T b, T c, T d) {
      t /= d;
      return c*t*t*t + b;
    }



    // Ord3 easing Out - decelerating to zero velocity

    template<typename T>
    static T easeOutOrd3(T t, T b, T c, T d) {
      t /= d;
      t--;
      return c*(t*t*t + static_cast<T>(1)) + b;
    }



    // Ord3 easing In/Out - acceleration until halfway, then deceleration

    template<typename T>
    static T easeInOutOrd3(T t, T b, T c, T d) {
      t /= d/static_cast<T>(2);
      if (t < static_cast<T>(1)) return c/static_cast<T>(2)*t*t*t + b;
      t -= static_cast<T>(2);
      return c/static_cast<T>(2)*(t*t*t + static_cast<T>(2)) + b;
    }


    // quartic easing In - accelerating from zero velocity

    template<typename T>
    static T easeInOrd4(T t, T b, T c, T d) {
      t /= d;
      return c*t*t*t*t + b;
    }



    // quartic easing Out - decelerating to zero velocity

    template<typename T>
    static T easeOutOrd4(T t, T b, T c, T d) {
      t /= d;
      t--;
      return -c * (t*t*t*t - static_cast<T>(1)) + b;
    }



    // quartic easing In/Out - acceleration until halfway, then deceleration

    template<typename T>
    static T easeInOutOrd4(T t, T b, T c, T d) {
      t /= d/static_cast<T>(2);
      if (t < static_cast<T>(1)) return c/static_cast<T>(2)*t*t*t*t + b;
      t -= static_cast<T>(2);
      return -c/static_cast<T>(2) * (t*t*t*t - static_cast<T>(2)) + b;
    }


    // quintic easing In - accelerating from zero velocity

    template<typename T>
    static T easeInOrd5(T t, T b, T c, T d) {
      t /= d;
      return c*t*t*t*t*t + b;
    }



    // quintic easing Out - decelerating to zero velocity

    template<typename T>
    static T easeOutOrd5(T t, T b, T c, T d) {
      t /= d;
      t--;
      return c*(t*t*t*t*t + static_cast<T>(1)) + b;
    }



    // quintic easing In/Out - acceleration until halfway, then deceleration

    template<typename T>
    static T easeInOutOrd5(T t, T b, T c, T d) {
      t /= d/static_cast<T>(2);
      if (t < static_cast<T>(1)) return c/static_cast<T>(2)*t*t*t*t*t + b;
      t -= static_cast<T>(2);
      return c/static_cast<T>(2)*(t*t*t*t*t + static_cast<T>(2)) + b;
    }


    // sinusoidal easing In - accelerating from zero velocity

    template<typename T>
    static T easeInSine(T t, T b, T c, T d) {
      return -c * std::cos(t/d * static_cast<T>(M_PI_2)) + c + b;
    }



    // sinusoidal easing Out - decelerating to zero velocity

    template<typename T>
    static T easeOutSine(T t, T b, T c, T d) {
      return c * std::sin(t/d * static_cast<T>(M_PI_2)) + b;
    }



    // sinusoidal easing In/Out - accelerating until halfway, then decelerating

    template<typename T>
    static T easeInOutSine(T t, T b, T c, T d) {
      return -c/static_cast<T>(2) * (std::cos(static_cast<T>(M_PI)*t/d) - static_cast<T>(1)) + b;
    }



    // exponential easing In - accelerating from zero velocity

    template<typename T>
    static T easeInExpo(T t, T b, T c, T d) {
      return c * std::pow( static_cast<T>(2), static_cast<T>(10) * (t/d - static_cast<T>(1)) ) + b;
    }

    template<typename T>
    static T proportional_value_derivative(T t_start_to_cur, T val_start, T val_range, T t_start_to_end) {

        // interpolation to compensate for "a small change in a small value has more effect than a small change in a big value"
        //
        // derivative of f should be proportional to f :
        //
        // f' = C f
        //
        // the solution has the form :
        //
        // f( x ) = - K exp( C x )
        //
        // initial / final conditions :
        //
        //      with :
        //
        //      t0 = time
        //      t1 = t0 + transition_duration
        //
        // f( t0 ) = cur
        // f( t1 ) = next
        //
        // give us a system of 2 equations with 2 unknowns leading to:
        //
        // C = - ln( f( t0 ) / f( t1 ) ) / ( t1 - t0 )
        // K = - f( t0 ) exp( - C t0 )
        //
        // so, with i := (t - t0) / (t1 - t0)
        //
        //      f(i) = f(t0) exp( -ln( f(t0)/f(t1) ) * i )
        // =>   f(i) = f(t0) exp(  ln( f(t1)/f(t0) ) * i )
        // =>   f(i) = f(t0) ( f(t1)/f(t0) ) ^ i )          (using a^(b*c) = a^b^c)
        //
        // which tranlates, with : i     -> t_start_to_cur / t_start_to_end
        //                         f(t0) -> val_start
        //                         f(t1) -> val_start + val_range
        // as:

        if( val_start == static_cast<T>(0) ) {
            return linearTween(t_start_to_cur, val_start, val_range, t_start_to_end);
        }
        auto val_end = val_start + val_range;
        if( val_end == static_cast<T>(0) ) {
            return linearTween(t_start_to_cur, val_start, val_range, t_start_to_end);
        }

        assert(same_sign_strict( val_end, val_start ));

        return val_start * std::pow(val_end/val_start, t_start_to_cur/t_start_to_end);
    }

    // exponential easing Out - decelerating to zero velocity

    template<typename T>
    static T easeOutExpo(T t, T b, T c, T d) {
      return c * ( -std::pow( static_cast<T>(2), -static_cast<T>(10) * t/d ) + static_cast<T>(1) ) + b;
    }



    // exponential easing In/Out - accelerating until halfway, then decelerating

    template<typename T>
    static T easeInOutExpo(T t, T b, T c, T d) {
      t /= d/static_cast<T>(2);
      if (t < static_cast<T>(1)) return c/static_cast<T>(2) * std::pow( static_cast<T>(2), static_cast<T>(10) * (t - static_cast<T>(1)) ) + b;
      t--;
      return c/static_cast<T>(2) * ( -std::pow( static_cast<T>(2), -static_cast<T>(10) * t) + static_cast<T>(2) ) + b;
    }


    // circular easing In - accelerating from zero velocity

    template<typename T>
    static T easeInCirc(T t, T b, T c, T d) {
      t /= d;
      return -c * (std::sqrt(static_cast<T>(1) - t*t) - static_cast<T>(1)) + b;
    }



    // circular easing Out - decelerating to zero velocity

    template<typename T>
    static T easeOutCirc(T t, T b, T c, T d) {
      t /= d;
      t--;
      return c * std::sqrt(static_cast<T>(1) - t*t) + b;
    }



    // circular easing In/Out - acceleration until halfway, then deceleration

    template<typename T>
    static T easeInOutCirc(T t, T b, T c, T d) {
      t /= d/static_cast<T>(2);
      if (t < static_cast<T>(1)) return -c/static_cast<T>(2) * (std::sqrt(static_cast<T>(1) - t*t) - static_cast<T>(1)) + b;
      t -= static_cast<T>(2);
      return c/static_cast<T>(2) * (std::sqrt(static_cast<T>(1) - t*t) + static_cast<T>(1)) + b;
    }

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
