
namespace imajuscule
{
    enum class FilterType
    {
        HIGH_PASS,
        LOW_PASS
    };
    
    template<FilterType KIND, int NDIMS, typename T>
    struct FilterState;
    
    template<int NDIMS, typename T>
    struct FilterState<FilterType::LOW_PASS, NDIMS, T> {
        using Tr = NumTraits<T>;
        void update(T alpha, T raw[NDIMS]) {
            for(int i=0; i<NDIMS; i++) {
                m_cur[i] = raw[i] * alpha + m_cur[i] * (Tr::one() - alpha);
            }
        }
        T m_cur[NDIMS]{};
    };

    template<int NDIMS, typename T>
    struct FilterState<FilterType::HIGH_PASS, NDIMS, T> {
        using Tr = NumTraits<T>;
        void update(T alpha, T raw[NDIMS]) {
            for(int i=0; i<NDIMS; i++) {
                m_cur[i] = alpha * (m_cur[i] + raw[i] - m_last[i]);
            }
            std::memcpy(m_last,raw,sizeof(m_last));
        }
        T m_cur[NDIMS]{};
        T m_last[NDIMS]{};
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

    template <class T, int NDIMS, FilterType KIND, bool ADAPTATIVE=false>
    class Filter
    {
        using Tr = NumTraits<T>;

        static constexpr auto kAccelerometerMinStep	= 0.02f;
        static constexpr auto kAccelerometerNoiseAttenuation = 3.0f;

    public:
        
        T square_magnitude(T rate, T freq) const {
            auto fcut = freq_from_cst(rate);
            A(fcut != 0);
            auto ratio = freq/fcut;
            return Tr::one() / get_inv_square_filter_magnitude<KIND>(ratio*ratio);
        }
        
        void initWithFreq(T rate, T cutOffFreq) {
            initWithAngleIncrement(Tr::two() * cutOffFreq / rate);
        }

        void initWithAngleIncrement(T inc) {
            auto im = inc * static_cast<T>(M_PI);
            
            if(KIND == FilterType::LOW_PASS) {
                FilterConstant = im / (im + Tr::one());
            }
            else if(KIND == FilterType::HIGH_PASS) {
                FilterConstant = Tr::one() / (Tr::one() + im);
            }
        }
        
        void feed(T raw[NDIMS]) {
            float alpha;
            if (ADAPTATIVE)
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
            else {
                alpha = FilterConstant;
            }
            state.update(alpha, raw);
        }
        
        const T * filtered() const {
            return state.m_cur;
        }
    private:
        FilterState<KIND, NDIMS, T> state;
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
        
        T freq_from_cst(T rate) const {
            if(KIND == FilterType::LOW_PASS) {
                return rate/((2.f * M_PI)*((1/FilterConstant) - 1));
            }
            else if(KIND == FilterType::HIGH_PASS) {
                return (rate/(2.f * M_PI)) * ((1 / FilterConstant) - 1);
            }
        }
    };
}
