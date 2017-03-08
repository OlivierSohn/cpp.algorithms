
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
    
    template <class T, int NDIMS, FilterType KIND, bool ADAPTATIVE=false>
    class Filter
    {
        using Tr = NumTraits<T>;

        static constexpr auto kAccelerometerMinStep	= 0.02f;
        static constexpr auto kAccelerometerNoiseAttenuation = 3.0f;

    public:
        void initWithSampleRate(T rate, T cutOffFreq) {
            // http://www.electronics-tutorials.ws/filter/filter_2.html : cutoff freq is where gain is -3db
            auto dt = Tr::one() / rate;
            auto RC = Tr::one() / (cutOffFreq * (2.f * M_PI));
            
            if(KIND == FilterType::LOW_PASS) {
                FilterConstant = dt / (dt + RC);
            }
            else if(KIND == FilterType::HIGH_PASS) {
                FilterConstant = RC / (dt + RC);
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
    };
}
