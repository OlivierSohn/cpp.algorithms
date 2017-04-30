/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule
{
    template < typename T > struct NumTraits; // not implemented

    template <class T>
    struct BaseNumTraits
    {
        using Type = T;
    };

    template <>
    struct NumTraits<float> : public BaseNumTraits<float> {
        using BaseNumTraits::Type;
        
        static constexpr Type zero() { return 0.f; }
        static constexpr Type half() { return .5f; }
        static constexpr Type one() { return 1.f; }
        static constexpr Type one_and_half() { return 1.5f; }
        static constexpr Type two() { return 2.f; }
        static constexpr Type three() { return 3.f; }
        
        static constexpr int discretize(Type fVal)
        {
            fVal += 0.5f;
            // cast float -> int rounds towards 0
            if (fVal < 0.f) {
                fVal -= 1.f;
            }
            return static_cast<int>(fVal);
        }
        
        // convert 3.00000 to 3 / 5. to 5 / 0.4000 to 0.4 / leave 3.00001 as is
        static std::string toSignificantString(const float val) {
            auto str = std::to_string(val);
            auto pos_point = str.find_first_of('.');
            if(pos_point != std::string::npos) {
                auto last_significant = str.find_last_not_of('0');
                if( last_significant != std::string::npos ) {
                    if(last_significant == pos_point) {
                        str.erase(pos_point);
                    }
                    else {
                        str.erase(last_significant + 1);
                    }
                }
            }
            return str;
        }
        
        static std::string to_string_with_precision(const float a_value, const int n)
        {
            std::ostringstream out;
            out << std::setprecision(n) << a_value;
            return out.str();
        }
    };
    
    template <>
    struct NumTraits<double> : public BaseNumTraits<double> {
        using BaseNumTraits::Type;
        
        static constexpr Type zero() { return 0.; }
        static constexpr Type half() { return .5; }
        static constexpr Type one() { return 1.; }
        static constexpr Type one_and_half() { return 1.5; }
        static constexpr Type two() { return 2.; }
        static constexpr Type three() { return 3.; }

        static constexpr int discretize(Type v)
        {
            v += 0.5;
            // cast -> int rounds towards 0
            if (v < 0.) {
                v -= 1.;
            }
            return static_cast<int>(v);
        }
    };
    
    template <>
    struct NumTraits<int> : public BaseNumTraits<int> {
        using BaseNumTraits::Type;
        
        static constexpr Type zero() { return 0; }
        static constexpr Type one() { return 1; }
        
        static constexpr Type discretize(Type i) { return i; }
    };
    
    template <>
    struct NumTraits<int16_t> : public BaseNumTraits<int16_t> {
        using BaseNumTraits::Type;
        
        static constexpr Type zero() { return 0; }
        static constexpr Type one() { return 1; }

        static constexpr Type discretize(Type i) { return i; }
    };
    template <>
    struct NumTraits<unsigned long> : public BaseNumTraits<unsigned long> {
        using BaseNumTraits::Type;
        
        static constexpr Type zero() { return 0; }
        static constexpr Type one() { return 1; }

        static constexpr Type discretize(Type i) { return i; }
    };
    
    template <>
    struct NumTraits<bool> : public BaseNumTraits<bool> {
        using BaseNumTraits::Type;
        
        static constexpr Type zero() { return false; }
        static constexpr Type one() { return true; }
    };

}
