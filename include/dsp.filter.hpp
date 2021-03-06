/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */
namespace imajuscule::audio {

    enum class FilterType
    {
        HIGH_PASS,
        LOW_PASS
    };

    template<FilterType KIND, int NDIMS, typename T>
    struct OrderState;

    template<int NDIMS, typename T>
    struct OrderState<FilterType::LOW_PASS, NDIMS, T> {
        using Tr = NumTraits<T>;
        OrderState() {
            reset();
        }

        void reset() {
            m_cur.fill({});
        }

        void setInitialValue(T val) {
          m_cur.fill(val);
        }

        void update(T const alpha, T const raw[NDIMS]) {
            for(int i=0; i<NDIMS; i++) {
                m_cur[i] = raw[i] * alpha + m_cur[i] * (Tr::one() - alpha);
            }
        }
        std::array<T, NDIMS> m_cur;
    };

    template<int NDIMS, typename T>
    struct OrderState<FilterType::HIGH_PASS, NDIMS, T> {
        using Tr = NumTraits<T>;
        OrderState() {
            reset();
        }

        void reset() {
            m_cur.fill({});
            m_last.fill({});
        }

        void setInitialValue(T val) {
          m_cur.fill(val);
          m_last.fill(val);
        }

        void update(T const alpha, T const raw[NDIMS]) {
            for(int i=0; i<NDIMS; i++) {
                m_cur[i] = alpha * (m_cur[i] + raw[i] - m_last[i]);
            }
            std::copy(raw, raw + NDIMS, std::begin(m_last));
        }
        std::array<T, NDIMS> m_cur, m_last;
    };

    template<FilterType KIND, typename T>
    T get_inv_square_filter_magnitude(T square_ratio) {
        if(KIND == FilterType::LOW_PASS) {
            return 1 + square_ratio;
        }
        else if(KIND == FilterType::HIGH_PASS) {
            return 1 + 1/square_ratio;
        }
    }

    template<typename orderState, int NORDER>
    struct FixedOrderFilter {
        static_assert(NORDER >= 1);

        constexpr int getOrder() const {
            return NORDER;
        }

        std::array<orderState, NORDER> m_orders;
    };

    template<typename orderState>
    struct VarOrderFilter_ {
        void setOrder(int n) {
            assert(n > 0);
            m_orders.resize(n);
        }

        int getOrder() const {
            return static_cast<int>(m_orders.size());
        }

        std::vector<orderState> m_orders;
    };

    template <class T, int NDIMS, FilterType KIND, typename Base>
    class FilterBase : public Base
    {
        using Base::m_orders;
        using Base::getOrder;

        using Tr = NumTraits<T>;

        //static constexpr auto kAccelerometerMinStep	= 0.02f;
        //static constexpr auto kAccelerometerNoiseAttenuation = 3.0f;
    public:

        // not tested
        T square_magnitude(T rate, T freq) const {
            auto fcut = freq_from_cst(rate);
            Assert(fcut != 0);
            auto ratio = freq/fcut;
            return Tr::one() / get_inv_square_filter_magnitude<KIND>(ratio * ratio);
        }

        void forgetPastSignals() {
            for(auto & o : m_orders) {
                o.reset();
            }
        }

        void setInitialValue(T val) {
          for(auto & o : m_orders) {
            o.setInitialValue(val);
          }
        }

        void initWithFreq(T rate, T cutOffFreq) {
            initWithAngleIncrement(Tr::two() * cutOffFreq / rate);
        }

        void initWithAngleIncrement(T inc) {
            auto im = inc * static_cast<T>(M_PI);
            auto denom = im + Tr::one();
            assert(denom);

            if constexpr (KIND == FilterType::LOW_PASS) {
                FilterConstant = im / denom;
            }
            else if constexpr (KIND == FilterType::HIGH_PASS) {
                FilterConstant = Tr::one() / denom;
            }
            assert(FilterConstant == FilterConstant);
        }

        void feed(T const raw[NDIMS]) {
            T alpha;
            // this template parameter was removed as it was not used
/*            if (ADAPTATIVE)
            {
                T d = Clamp(fabs(Norm(state.m_cur) - Norm(raw)) / (T)kAccelerometerMinStep - Tr::one(),
                            Tr::zero(),
                            Tr::one());

                if(KIND == FilterType::LOW_PASS)
                {
                    alpha = (Tr::one() - d) * FilterConstant / (T)kAccelerometerNoiseAttenuation + d * FilterConstant;
                }
                else if(KIND == FilterType::HIGH_PASS)
                {
                    alpha = d * FilterConstant / (T)kAccelerometerNoiseAttenuation + (Tr::one() - d) * FilterConstant;
                }
                else {
                    assert(0);
                }
            }
            else */
            {
                alpha = FilterConstant;
            }

            {
                auto array = raw;
                for(auto & o : m_orders) {
                    o.update(alpha, array);
                    array = &o.m_cur[0];
                }
            }
        }

        const T * filtered(int order = 0) const {
            int const o = order ? order : (getOrder()-1);
            return &m_orders[o].m_cur[0];
        }

      T getConstant() const {
        return FilterConstant;
      }
    private:
        T FilterConstant;

        static T Norm(T v[NDIMS])
        {
            auto val = Tr::zero();
            for(int i=0;i<NDIMS;i++)
            {
                val += v[i]*v[i];
            }
            return sqrt(val);
        }

        static T Clamp(T v, T min, T max)
        {
            if(v > max)
                return max;
            else if(v < min)
                return min;
            else
                return v;
        }

        // not tested
        T freq_from_cst(T rate) const {
            if(KIND == FilterType::LOW_PASS) {
                return rate/((2.f * M_PI)*((1/FilterConstant) - 1));
            }
            else if(KIND == FilterType::HIGH_PASS) {
                return (rate/(2.f * M_PI)) * ((1 / FilterConstant) - 1);
            }
        }
    };

    template <class T, int NDIMS, FilterType KIND, int NORDER>
    struct FilterSelector {
        using type = FilterBase<T, NDIMS, KIND, FixedOrderFilter<OrderState<KIND, NDIMS, T>, NORDER>>;
    };

    template <class T, int NDIMS, FilterType KIND>
    struct FilterSelector<T, NDIMS, KIND, 0> {
        using type = FilterBase<T, NDIMS, KIND, VarOrderFilter_<OrderState<KIND, NDIMS, T>>>;
    };

    static constexpr auto VariableOrder = 0;

    template <class T, int NDIMS, FilterType KIND, int NORDER=1>
    using Filter = typename FilterSelector<T, NDIMS, KIND, NORDER>::type;
}
