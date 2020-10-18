
namespace imajuscule
{
    // Inverse of sRGB "gamma" function. (approx 2.2)
    static inline float inv_gam_sRGB(float c) {
        if ( c <= 0.04045f )
            return c/12.92f;
        else
            return pow(((c+0.055f)/(1.055f)),2.4f);
    }
    
    // sRGB "gamma" function (approx 2.2)
    static inline float gam_sRGB(float v) {
        if(v<=0.0031308f)
            return v * 12.92f;
        else
            return 1.055f*pow(v,1.0f/2.4f)-0.055f;
    }
    
    static inline unsigned char linearToSRGB(float c) {
        auto res = inv_gam_sRGB(c);
        assert(res <= 1.f);
        assert(res >= 0.f);
        return (unsigned char)(res*255.f + 0.5f);
    }
    
    static inline float SRGBToLinear(unsigned char ref) {
        float c = ref/255.f;
        auto res = gam_sRGB(c);
        assert(res <= 1.f);
        assert(res >= 0.f);
        return res;
    }
    
    static inline std::array<float, 3> linearToHSV(std::array<float, 3> const & ref) {
        // using http://www.rapidtables.com/convert/color/rgb-to-hsv.htm
        
        float r = ref[0];
        float g = ref[1];
        float b = ref[2];

        int maxComponent = 0;
        float cmax = r; // note that v == cmax
        if(g > cmax) {
            maxComponent = 1;
            cmax = g;
        }
        if(b > cmax) {
            maxComponent = 2;
            cmax = b;
        }
        float cmin = std::min(r,std::min(g,b));
        auto delta = cmax-cmin;

        if (delta < 0.00001) {
            // grey case
            return {{0,0,cmax}};
        }
        assert(cmax > 0);
        auto s = delta / cmax;
        float h;
        switch(maxComponent) {
            case 0:
                h = (g-b) / delta;
                break;
            case 1:
                h = 2.f + (b-r)/delta;
                break;
          default:
            case 2:
                h = 4.f + (r-g)/delta;
                break;
        }
        h /= 6.f; // normalize
        return {{h,s,cmax}};
    }
    
    template<typename T>
    struct LazyPlainColor {
        using value_type = T;
        
        LazyPlainColor(std::array<value_type, 3> && colors) :
        color(std::move(colors))
        , uptodate(true) {
        }
        
        LazyPlainColor() {
            uptodate = false;
        }
        
        bool isValid() const { return uptodate; }
        void invalidate() { uptodate = false; }
        
        template<typename F>
        void set(std::array<unsigned char,3> const & ref, F conversion) {
            for(int i=0; i<static_cast<int>(color.size()); ++i) {
                color[i] = conversion(ref[i]);
            }
            uptodate = true;
        }
        
        template<typename F>
        void setDiscrete(std::array<float,3> const & ref, F conversion) {
            for(int i=0; i<static_cast<int>(color.size()); ++i) {
                color[i] = conversion(ref[i]);
            }
            uptodate = true;
        }
        
        template<typename F, typename U>
        void set3(std::array<U,3> const & ref, F conversion) {
            color = conversion(ref);
            uptodate = true;
        }
        
        auto const & get() const {
            assert(uptodate);
            return color;
        }
        
        bool operator == (const LazyPlainColor & other) const {
            if(!uptodate && !other.uptodate) {
                return true;
            }
            if(!uptodate || !other.uptodate) {
                return false;
            }
            return color == other.color;
        }
    private:
        
        std::array<value_type,3> color;
        bool uptodate : 1;
    };
    
    struct Color8 {
        Color8() = default;
        
        Color8(std::array<unsigned char, 3> && colors) :
        srgb(std::move(colors)) {
        }
        
        Color8(const unsigned char arr[3]) :
        Color8(std::array<unsigned char, 3>{{
            arr[0],
            arr[1],
            arr[2]
        }}) {
        }
        
        void reset() {
            srgb.invalidate();
            linear.invalidate();
            hsv.invalidate();
        }
        
        // note this is NOT "color equality", it is "C++ object content equality".
        bool operator == (const Color8 & other) const {
            return (srgb == other.srgb) && (linear == other.linear);
        }
        
        void setInLinear(std::array<float, 3> && a) {
            srgb.invalidate();
            hsv.invalidate();
            linear = {std::move(a)};
        }
        void setInSRGB(std::array<unsigned char, 3> && a) {
            linear.invalidate();
            hsv.invalidate();
            srgb = {std::move(a)};
        }
        std::array<unsigned char,3> const & inSRGB() {
            if(!srgb.isValid()) {
                srgb.setDiscrete(linear.get(), linearToSRGB);
            }
            return srgb.get();
        }
        std::array<float,3> const & inLinearRGB() {
            if(!linear.isValid()) {
                linear.set(srgb.get(), SRGBToLinear);
            }
            return linear.get();
        }
        std::array<float,3> const & inHSV() {
            if(!hsv.isValid()) {
                hsv.set3(inLinearRGB(), linearToHSV);
            }
            return hsv.get();
        }
    private:
        // this is the only colorspace using discretized colors because it is given as input in discretized mode
        LazyPlainColor<unsigned char> srgb;
        LazyPlainColor<float> linear;
        LazyPlainColor<float> hsv;
    };
    
    void adjustColorBrightness(Color8 &, float brightness, float effectRatio);

    float colorBrightness(Color8 &);
    float squaredEuclidianDistance(Color8 &, Color8 &);
    float stackOverflowSquaredHSVDistance(Color8 &, Color8 &);
    
}
