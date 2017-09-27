
namespace imajuscule {
    
    void alphaCompositeColor(Color8 & color, float overlayed, float overlayAlpha) {
        auto const & c = color.inSRGB();
        std::array<unsigned char, 3> srgb;
        
        for(int i=0; i<3; ++i) {
            float v = (c[i]/255.f) * (1.f - overlayAlpha) + overlayed * overlayAlpha;
            srgb[i] = (unsigned char)(255.f * v);
        }
        
        color.setInSRGB(std::move(srgb));
    }
    
    
    float colorBrightness(Color8 & color) {
        // sRGB luminance(Y) values
        constexpr auto rY = 0.212655f;
        constexpr auto gY = 0.715158f;
        constexpr auto bY = 0.072187f;
        
        // perceived
        /*
        constexpr auto rY = 0.299f;
        constexpr auto gY = 0.587f;
        constexpr auto bY = 0.114f;
         */
        /*
        constexpr auto rY = 0.333f;
        constexpr auto gY = 0.333f;
        constexpr auto bY = 0.333f;
        */
        auto const & c = color.inSRGB();
        
        return
        rY * (c[0]/255.f) +
        gY * (c[1]/255.f) +
        bY * (c[2]/255.f);
    }
    
    void adjustColorBrightness(Color8 & color, float targetBrightness, float effectRatio) {
        using namespace std;
        auto const & c = color.inSRGB();
        unsigned char M = std::max(c[0], std::max(c[1], c[2]));
        if(!M) {
            return;
        }
        
        auto brightness = colorBrightness(color);
        assert(brightness);
        
        targetBrightness = brightness + effectRatio * (targetBrightness-brightness);
        
        if(brightness == targetBrightness) {
            return;
        }
        
        // mixed approach : scale until one component saturates, then alpha compose with white
        auto scaleRatioMax = 255.5f / (float)M;
        auto generalRatio = targetBrightness / brightness;
        bool needAlphaCompositing = generalRatio > scaleRatioMax;
        auto scaleRatio = std::min(scaleRatioMax, generalRatio);
        
        std::array<unsigned char, 3> srgbCandidate;
        
        for(int i=0; i<3; ++i) {
            // it is easy to match brightness this way
            // but it might look better to do alpha compositing with white instead
            float res = c[i] * scaleRatio;
            assert(res < 255.8f);
            srgbCandidate[i] = res;
        }
        color.setInSRGB(std::move(srgbCandidate));
        
        if(!needAlphaCompositing) {
            return;
        }
        // update brightness
        brightness = colorBrightness(color);
        
        float alpha;
        float overlay;
        if(brightness < targetBrightness) {
            // white is 1
            overlay = 1.f;
            
            // c is brightness
            // alpha compositing: targetBrightness = (1-alpha) * brightness + alpha * 1
            // alpha * (1 - brightness) = targetBrightness - brightness
            alpha = 1.f - ((1.f - targetBrightness) / (1.f-brightness));
        }
        else {
            // black is 0
            overlay = 0.f;
            
            // c is brightness
            // alpha compositing: targetBrightness = (1-alpha) * brightness
            // alpha * (- brightness) = targetBrightness - brightness
            alpha = 1.f - (targetBrightness / brightness);
            
            // note that this branch is equivalent to Scaling method
        }
        assert(alpha <= 1.f);
        assert(alpha >= 0.f);
        alphaCompositeColor(color, overlay, alpha);
        //cout << "brightness: " << setprecision(2) << brightness << " " << colorBrightness(color) << " " << targetBrightness << endl;
    }
 
    float stackOverflowSquaredHSVDistance(Color8 & color1, Color8 & color2) {
        // cf Sean Gerrish answer https://stackoverflow.com/questions/35113979/calculate-distance-between-colors-in-hsv-space
        auto const c1 = color1.inHSV();
        auto const c2 = color2.inHSV();
        auto const h1 = c1[0];
        auto const h2 = c2[0];
        auto const s1 = c1[1];
        auto const s2 = c2[1];
        auto const v1 = c1[2];
        auto const v2 = c2[2];

        auto const tmp1 = std::sin(2.f*(float)M_PI*h1)*s1*v1 - std::sin(2.f*(float)M_PI*h2)*s2*v2;
        auto const tmp2 = std::cos(2.f*(float)M_PI*h1)*s1*v1 - std::cos(2.f*(float)M_PI*h2)*s2*v2;
        auto const tmp3 = v1-v2;
        return tmp1*tmp1
        + tmp2*tmp2
        + tmp3*tmp3;
    }

    float squaredEuclidianDistance(Color8 & color1, Color8 & color2) {
        // https://en.wikipedia.org/wiki/Color_difference
        auto c1 = color1.inSRGB();
        auto c2 = color2.inSRGB();
        
        auto r = (c1[0] + c2[0]) / (255.f * 2.f);
        auto dR = c1[0] - c2[0];
        auto dG = c1[1] - c2[1];
        auto dB = c1[2] - c2[2];
        
        return (dR*dR*(2 + r) + 4*dG*dG + dB*dB*(3 - r)) / (255.f * 255.f);
    }
}
