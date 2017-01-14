

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
        static constexpr Type zero() { return 0.f; }
        static constexpr Type half() { return .5f; }
        static constexpr Type one() { return 1.f; }
        static constexpr Type two() { return 2.f; }
        static constexpr Type three() { return 3.f; }
    };
    
    template <>
    struct NumTraits<double> : public BaseNumTraits<double> {
        static constexpr Type zero() { return 0.; }
        static constexpr Type half() { return .5; }
        static constexpr Type one() { return 1.; }
        static constexpr Type two() { return 2.; }
        static constexpr Type three() { return 3.; }
    };
    
    template <>
    struct NumTraits<int> : public BaseNumTraits<int> {
        static constexpr Type zero() { return 0; }
        static constexpr Type one() { return 1; }
        
    };
    
    template <>
    struct NumTraits<int16_t> : public BaseNumTraits<int16_t> {
        static constexpr Type zero() { return 0; }
        static constexpr Type one() { return 1; }
    };
    template <>
    struct NumTraits<unsigned long> : public BaseNumTraits<unsigned long> {
        static constexpr Type zero() { return 0; }
        static constexpr Type one() { return 1; }
    };
    
    template <>
    struct NumTraits<bool> : public BaseNumTraits<bool> {
        static constexpr Type zero() { return false; }
        static constexpr Type one() { return true; }
    };

}
