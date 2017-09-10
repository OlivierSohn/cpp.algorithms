
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
    
    static inline unsigned char linearToSRGB(unsigned char ref) {
        float c = ref/255.f;
        auto res = inv_gam_sRGB(c);
        assert(res <= 1.f);
        assert(res >= 0.f);
        return (unsigned char)(res*255 + 0.5f);
    }
    
    static inline unsigned char SRGBToLinear(unsigned char ref) {
        float c = ref/255.f;
        auto res = gam_sRGB(c);
        assert(res <= 1.f);
        assert(res >= 0.f);
        return (unsigned char)(res*255 + 0.5f);
    }
    
    struct LazyPlainColor {
        LazyPlainColor(std::array<unsigned char, 3> && colors) :
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
            for(int i=0; i<color.size(); ++i) {
                color[i] = conversion(ref[i]);
            }
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
        
        std::array<unsigned char,3> color;
        bool uptodate : 1;
    };
    
    struct Color8 {
        Color8() = default;
        
        Color8(std::array<unsigned char, 3> && colors) : srgb(std::move(colors)) {
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
        }
        
        // note this is NOT "color equality", it is "C++ object content equality".
        bool operator == (const Color8 & other) const {
            return (srgb == other.srgb) && (linear == other.linear);
        }
        
        void setInLinear(std::array<unsigned char, 3> && a) {
            srgb.invalidate();
            linear = {std::move(a)};
        }
        void setInSRGB(std::array<unsigned char, 3> && a) {
            linear.invalidate();
            srgb = {std::move(a)};
        }
        std::array<unsigned char,3> const & inSRGB() {
            if(!srgb.isValid()) {
                srgb.set(linear.get(), linearToSRGB);
            }
            return srgb.get();
        }
        std::array<unsigned char,3> const & inLinearRGB() {
            if(!linear.isValid()) {
                linear.set(srgb.get(), SRGBToLinear);
            }
            return linear.get();
        }
    private:
        LazyPlainColor srgb, linear;
    };
    
    
    enum class BrightnessAdjustment {
        AlphaCompositing,
        Scaling
    };
    
    void adjustColorBrightness(BrightnessAdjustment, Color8 &, float brightness);

    float colorBrightnessLinear(Color8 &);
    float squaredEuclidianDistance(Color8 &, Color8 &);
    
}
