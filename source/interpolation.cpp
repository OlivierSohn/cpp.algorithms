/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

using namespace imajuscule;

namespace imajuscule {

    // we add 1 for every group, because a value of 0.f means the first value of the group
    const float itp::NO_EASE = safe_cast<float>(INTERPOLATION_LOWER_BOUND + 1);
    const float itp::EASE_IN = safe_cast<float>(EASE_IN_first_value + 1);
    const float itp::EASE_OUT = safe_cast<float>(EASE_OUT_first_value + 1);
    const float itp::EASE_IO = safe_cast<float>(EASE_INOUT_first_value + 1);

    // values for NO_EASE:
    const float itp::P1 = 0.f;
    const float itp::PVD = 1.f;

    // values for EASE_IN/OUT/IO:
    const float itp::P2 = 0.f;
    const float itp::P3 = 1.f;
    const float itp::P4 = 2.f;
    const float itp::P5 = 3.f;
    const float itp::SINE = 4.f;
    const float itp::EXPO = 5.f;
    const float itp::CIRC = 6.f;
}
enumTraversal itp::gInterpolationTraversal(
    itp::INTERPOLATION_LOWER_BOUND,
    itp::INTERPOLATION_UPPER_BOUND,
    [](int val)->const char* {
    return itp::interpolationInfo(val);
},
[](int val)->bool {
    return itp::intIsReal(val);
}
);
const enumTraversal & itp::interpolation_traversal()
{
    return gInterpolationTraversal;
}

itp::interpolation itp::intToInterpolation(int val, bool & bReal)
{
    auto iVal = safe_cast<itp::interpolation>(val);
    bReal = intIsReal(val);
    return iVal;
}

const char * itp::interpolationInfo(int val)
{
    itp::interpolation itpVal = (itp::interpolation) val;

    switch (itpVal)
    {
    case LINEAR:
        return "Linear";
        break;
    case PROPORTIONAL_VALUE_DERIVATIVE:
        return "Proportional Value Derivative";
    case EASE_IN_QUAD:
        return "Ease In Quad";
        break;
    case EASE_IN_CUBIC:
        return "Ease In Ord3";
        break;
    case EASE_IN_QUART:
        return "Ease In Quart";
        break;
    case EASE_IN_QUINT:
        return "Ease In Quint";
        break;
    case EASE_IN_SINE:
        return "Ease In Sine";
        break;
    case EASE_IN_EXPO:
        return "Ease In Exp";
        break;
    case EASE_IN_CIRC:
        return "Ease In Circ";
        break;
    case EASE_OUT_QUAD:
        return "Ease Out Quad";
        break;
    case EASE_OUT_CUBIC:
        return "Ease Out Ord3";
        break;
    case EASE_OUT_QUART:
        return "Ease Out Quart";
        break;
    case EASE_OUT_QUINT:
        return "Ease Out Quint";
        break;
    case EASE_OUT_SINE:
        return "Ease Out Sine";
        break;
    case EASE_OUT_EXPO:
        return "Ease Out Exp";
        break;
    case EASE_OUT_CIRC:
        return "Ease Out Circ";
        break;
    case EASE_INOUT_QUAD:
        return "Ease InOut Quad";
        break;
    case EASE_INOUT_CUBIC:
        return "Ease InOut Ord3";
        break;
    case EASE_INOUT_QUART:
        return "Ease InOut Quart";
        break;
    case EASE_INOUT_QUINT:
        return "Ease InOut Quint";
        break;
    case EASE_INOUT_SINE:
        return "Ease InOut Sine";
        break;
    case EASE_INOUT_EXPO:
        return "Ease InOut Exp";
        break;
    case EASE_INOUT_CIRC:
        return "Ease InOut Circ";
        break;
    default:
        return "missing description";
        break;
    }
}

float itp::interpolate( interpolation type, float t, float b, float c, float d)
{
    switch(type)
    {
        default:
            assert(0);
        case LINEAR:
            return linearTween(t,b,c,d);
            break;
        case PROPORTIONAL_VALUE_DERIVATIVE:
            return proportional_value_derivative(t,b,c,d);
            break;
        case EASE_IN_QUAD:
            return easeInOrd2(t,b,c,d);
            break;
        case EASE_IN_CUBIC:
            return easeInOrd3(t,b,c,d);
            break;
        case EASE_IN_QUART:
            return easeInOrd4(t,b,c,d);
            break;
        case EASE_IN_QUINT:
            return easeInOrd5(t,b,c,d);
            break;
        case EASE_IN_SINE:
            return easeInSine(t,b,c,d);
            break;
        case EASE_IN_EXPO:
            return easeInExpo(t,b,c,d);
            break;
        case EASE_IN_CIRC:
            return easeInCirc(t,b,c,d);
            break;
        case EASE_OUT_QUAD:
            return easeOutOrd2(t,b,c,d);
            break;
        case EASE_OUT_CUBIC:
            return easeOutOrd3(t,b,c,d);
            break;
        case EASE_OUT_QUART:
            return easeOutOrd4(t,b,c,d);
            break;
        case EASE_OUT_QUINT:
            return easeOutOrd5(t,b,c,d);
            break;
        case EASE_OUT_SINE:
            return easeOutSine(t,b,c,d);
            break;
        case EASE_OUT_EXPO:
            return easeOutExpo(t,b,c,d);
            break;
        case EASE_OUT_CIRC:
            return easeOutCirc(t,b,c,d);
            break;
        case EASE_INOUT_QUAD:
            return easeInOutOrd2(t,b,c,d);
            break;
        case EASE_INOUT_CUBIC:
            return easeInOutOrd3(t,b,c,d);
            break;
        case EASE_INOUT_QUART:
            return easeInOutOrd4(t,b,c,d);
            break;
        case EASE_INOUT_QUINT:
            return easeInOutOrd5(t,b,c,d);
            break;
        case EASE_INOUT_SINE:
            return easeInOutSine(t,b,c,d);
            break;
        case EASE_INOUT_EXPO:
            return easeInOutExpo(t,b,c,d);
            break;
        case EASE_INOUT_CIRC:
            return easeInOutCirc(t,b,c,d);
            break;
    }
    assert(0);
    return 0.f;
}
//simple linear tweening - no easing, no acceleration



float itp::linearTween(
    float t, // t_cur - t_start
    float b, // val_start
    float c, // val_end - val_start
    float d  // t_end - t_start
    )
{
	return c*t/d + b;
}


// quadratic easing In - accelerating from zero velocity


float itp::easeInOrd2(float t, float b, float c, float d) {
	t /= d;
	return c*t*t + b;
}


// quadratic easing Out - decelerating to zero velocity


float itp::easeOutOrd2(float t, float b, float c, float d) {
	t /= d;
	return -c * t*(t-2.0f) + b;
}



// quadratic easing In/Out - acceleration until halfway, then deceleration


float itp::easeInOutOrd2(float t, float b, float c, float d) {
	t /= d/2.0f;
	if (t < 1.0f) return c/2.0f*t*t + b;
	t--;
	return -c/2.0f * (t*(t-2.0f) - 1.0f) + b;
}


// Ord3 easing In - accelerating from zero velocity


float itp::easeInOrd3(float t, float b, float c, float d) {
	t /= d;
	return c*t*t*t + b;
}



// Ord3 easing Out - decelerating to zero velocity


float itp::easeOutOrd3(float t, float b, float c, float d) {
	t /= d;
	t--;
	return c*(t*t*t + 1.0f) + b;
}



// Ord3 easing In/Out - acceleration until halfway, then deceleration


float itp::easeInOutOrd3(float t, float b, float c, float d) {
	t /= d/2.0f;
	if (t < 1.0f) return c/2.0f*t*t*t + b;
	t -= 2.0f;
	return c/2.0f*(t*t*t + 2.0f) + b;
}


// quartic easing In - accelerating from zero velocity


float itp::easeInOrd4(float t, float b, float c, float d) {
	t /= d;
	return c*t*t*t*t + b;
}



// quartic easing Out - decelerating to zero velocity


float itp::easeOutOrd4(float t, float b, float c, float d) {
	t /= d;
	t--;
	return -c * (t*t*t*t - 1.0f) + b;
}



// quartic easing In/Out - acceleration until halfway, then deceleration


float itp::easeInOutOrd4(float t, float b, float c, float d) {
	t /= d/2.0f;
	if (t < 1.0f) return c/2.0f*t*t*t*t + b;
	t -= 2.0f;
	return -c/2.0f * (t*t*t*t - 2.0f) + b;
}


// quintic easing In - accelerating from zero velocity


float itp::easeInOrd5(float t, float b, float c, float d) {
	t /= d;
	return c*t*t*t*t*t + b;
}



// quintic easing Out - decelerating to zero velocity


float itp::easeOutOrd5(float t, float b, float c, float d) {
	t /= d;
	t--;
	return c*(t*t*t*t*t + 1.0f) + b;
}



// quintic easing In/Out - acceleration until halfway, then deceleration


float itp::easeInOutOrd5(float t, float b, float c, float d) {
	t /= d/2.0f;
	if (t < 1.0f) return c/2.0f*t*t*t*t*t + b;
	t -= 2.0f;
	return c/2.0f*(t*t*t*t*t + 2.0f) + b;
}


// sinusoidal easing In - accelerating from zero velocity


float itp::easeInSine(float t, float b, float c, float d) {
	return -c * cosf(t/d * ((float)M_PI_2)) + c + b;
}



// sinusoidal easing Out - decelerating to zero velocity


float itp::easeOutSine(float t, float b, float c, float d) {
	return c * sinf(t/d * ((float)M_PI_2)) + b;
}



// sinusoidal easing In/Out - accelerating until halfway, then decelerating


float itp::easeInOutSine(float t, float b, float c, float d) {
	return -c/2.0f * (cosf((float)M_PI*t/d) - 1.0f) + b;
}



// exponential easing In - accelerating from zero velocity


float itp::easeInExpo(float t, float b, float c, float d) {
	return c * powf( 2.0f, 10.0f * (t/d - 1.0f) ) + b;
}

float itp::proportional_value_derivative(float t_start_to_cur, float val_start, float val_range, float t_start_to_end) {

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

    if( val_start == 0.f ) {
        return linearTween(t_start_to_cur, val_start, val_range, t_start_to_end);
    }
    auto val_end = val_start + val_range;
    if( val_end == 0.f ) {
        return linearTween(t_start_to_cur, val_start, val_range, t_start_to_end);
    }

    assert(same_sign_strict( val_end, val_start ));

    return val_start * powf(val_end/val_start, t_start_to_cur/t_start_to_end);
}

// exponential easing Out - decelerating to zero velocity


float itp::easeOutExpo(float t, float b, float c, float d) {
	return c * ( -powf( 2.0f, -10.0f * t/d ) + 1.0f ) + b;
}



// exponential easing In/Out - accelerating until halfway, then decelerating


float itp::easeInOutExpo(float t, float b, float c, float d) {
	t /= d/2.0f;
	if (t < 1.0f) return c/2.0f * powf( 2.0f, 10.0f * (t - 1.0f) ) + b;
	t--;
	return c/2.0f * ( -powf( 2.0f, -10.0f * t) + 2.0f ) + b;
}


// circular easing In - accelerating from zero velocity


float itp::easeInCirc(float t, float b, float c, float d) {
	t /= d;
	return -c * (sqrtf(1.0f - t*t) - 1.0f) + b;
}



// circular easing Out - decelerating to zero velocity


float itp::easeOutCirc(float t, float b, float c, float d) {
	t /= d;
	t--;
	return c * sqrtf(1.0f - t*t) + b;
}



// circular easing In/Out - acceleration until halfway, then deceleration


float itp::easeInOutCirc(float t, float b, float c, float d) {
	t /= d/2.0f;
	if (t < 1.0f) return -c/2.0f * (sqrtf(1.0f - t*t) - 1.0f) + b;
	t -= 2.0f;
	return c/2.0f * (sqrtf(1.0f - t*t) + 1.0f) + b;
}
