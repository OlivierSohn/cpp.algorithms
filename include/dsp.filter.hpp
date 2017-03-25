
namespace imajuscule
{
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
        
        void update(T const alpha, T raw[NDIMS]) {
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
        
        void update(T const alpha, T raw[NDIMS]) {
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
        static_assert(NORDER >= 1, "");

        int getOrder() const {
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
        
        static constexpr auto kAccelerometerMinStep	= 0.02f;
        static constexpr auto kAccelerometerNoiseAttenuation = 3.0f;
    public:
        
        // not tested
        T square_magnitude(T rate, T freq) const {
            auto fcut = freq_from_cst(rate);
            A(fcut != 0);
            auto ratio = freq/fcut;
            return Tr::one() / get_inv_square_filter_magnitude<KIND>(ratio * ratio);
        }
        
        void forgetPastSignals() {
            for(auto & o : m_orders) {
                o.reset();
            }
        }
        
        void initWithFreq(T rate, T cutOffFreq) {
            initWithAngleIncrement(Tr::two() * cutOffFreq / rate);
        }

        void initWithAngleIncrement(T inc) {
            auto im = inc * static_cast<T>(M_PI);
            auto denom = im + Tr::one();
            assert(denom);
                                      
            if(KIND == FilterType::LOW_PASS) {
                FilterConstant = im / denom;
            }
            else if(KIND == FilterType::HIGH_PASS) {
                FilterConstant = Tr::one() / denom;
            }
            assert(FilterConstant == FilterConstant);
        }
        
        void feed(T raw[NDIMS]) {
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
    
    
    template<typename T>
    struct FIRFilter {
        template<typename U>
        FIRFilter(std::vector<U> const & c) : past(c.size(), {}) {
            assert(c.size() == past.size());
            coefficients.reserve(c.size());
            coefficients.clear();
            for(auto coeff : c) {
                coefficients.push_back(coeff);
            }
        }
        
        auto size() const { return coefficients.size(); }
        
        void step(T val) {
            past.feed(val);
        }
        
        T get() const {
            auto res = T{};
            int index = 0;
            // when coefficients are symmetrical it doesn't matter
            // if we are traversing forward or backward
            past.for_each_bkwd([&res, &index, this](auto val) {
                res += val * coefficients[index];
                ++index;
            });
            return res;
        }

    private:
        std::vector<double> coefficients;
        cyclic<T> past;
    };
    
    
    static void plotMagnitude(fft::FFTVec<double> const & v) {
        std::vector<double> mags;
        mags.reserve(v.size());
        std::transform(v.begin(), v.end(),
                       std::back_inserter(mags),
                       [](auto v){return abs(v);});
        StringPlot plot(30,1024);
        plot.draw(mags);
        plot.log();
    }
    
    template<typename Container, typename T, typename F>
    auto sample_frequencies(Container & res, T nyquist_freq, F getFreq) {
    
        auto N = res.size();
        auto NumTaps = N;
        auto nyquist = NumTaps/2;
        
        double RadPerSample = -M_PI;
        if(0 == NumTaps % 2) {
            RadPerSample *= (N - 1.0)/N;
        }
        for(int i=0; i<=nyquist; ++i) {
            auto f = nyquist_freq * i / nyquist;
            T magnitude = getFreq(f);
            auto cplx = magnitude * polar(RadPerSample*i);
            res[i] = cplx;
            if(i && i != nyquist) {
                auto conji = N-i;
                res[conji] = cplx;
            }
        }

    }
    
    template<typename T, typename F>
    auto fir_coefficients_by_f_sampling(T nyquist_freq, F getFreq) {
        ScopedLog l("Compute", "FIR coeffs by freq. sampling");
        // according to http://iowahills.com/FIRFiltersByFreqSampling.html
        // with the same number of taps as of fft size (we could try less)
        
        using namespace imajuscule::fft;
        
        constexpr auto N = 4096;
        FFTVec<double> res, input;
        res.resize(N);
        input.resize(N);
        
        sample_frequencies(input, nyquist_freq, getFreq);

        compute_fft(N, input.begin(), res.begin());
        normalize_fft(N, res.begin(), res.end());
        apply_hann_window(res.begin(), res.end());
        
        //plotMagnitude(res);
        
        std::vector<T> v;
        v.reserve(res.size());
        std::transform(res.begin(), res.end(),
                       std::back_inserter(v),
                       [](auto val) { return val.real(); });
        return v;
    }
    
    
}
