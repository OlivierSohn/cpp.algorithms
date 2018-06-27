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
