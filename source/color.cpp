
namespace imajuscule {
    
    void alphaCompositeColor(Color8 & color, float overlayed, float overlayAlpha) {
        auto const & c = color.inLinearRGB();
        std::array<float, 3> linear;
        
        for(int i=0; i<3; ++i) {
            linear[i] = (c[i] * (1 - overlayAlpha) + overlayed * overlayAlpha);
        }
        
        color.setInLinear(std::move(linear));
    }
    
    
    float colorBrightnessLinear(Color8 & color) {
        // sRGB luminance(Y) values
        constexpr auto rY = 0.212655f;
        constexpr auto gY = 0.715158f;
        constexpr auto bY = 0.072187f;
        
        auto const & c = color.inLinearRGB();
        
        return
        rY * c[0] +
        gY * c[1] +
        bY * c[2];
    }
    
    void adjustColorBrightness(BrightnessAdjustment adj, Color8 & color, float targetLinearBrightness, float effectRatio) {
        using namespace std;
        auto linearBrightness = colorBrightnessLinear(color);
        if(!linearBrightness) {
            return;
        }
        
        targetLinearBrightness = linearBrightness + effectRatio * (targetLinearBrightness-linearBrightness);
        
        if(linearBrightness == targetLinearBrightness) {
            return;
        }
        
        if(adj == BrightnessAdjustment::AlphaCompositing) {
            float alpha;
            float overlay;
            if(linearBrightness < targetLinearBrightness) {
                // white is 1
                overlay = 1.f;
                
                // c is brightness
                // alpha compositing: targetBrightness = (1-alpha) * brightness + alpha * 1
                // alpha * (1 - brightness) = targetBrightness - brightness
                alpha = 1.f - ((1.f - targetLinearBrightness) / (1.f-linearBrightness));
            }
            else {
                // black is 0
                overlay = 0.f;
                
                // c is brightness
                // alpha compositing: targetBrightness = (1-alpha) * brightness
                // alpha * (- brightness) = targetBrightness - brightness
                alpha = 1.f - (targetLinearBrightness / linearBrightness);
                
                // note that this branch is equivalent to Scaling method
            }
            assert(alpha <= 1.f);
            assert(alpha >= 0.f);
            alphaCompositeColor(color, overlay, alpha);
        }
        else {
            auto ratio = targetLinearBrightness / linearBrightness;
            
            std::array<float, 3> linearCandidate;
            auto const & c = color.inLinearRGB();
            
            for(int i=0; i<3; ++i) {
                // it is easy to match brightness this way
                // but it might look better to do alpha compositing with white instead
                auto res = c[i] * ratio;
                if(res > 1) {
                    cout << "color brightness adjustment : saturation" << endl;
                    res = 1;
                }
                linearCandidate[i] = res;
            }
            color.setInLinear(std::move(linearCandidate));
        }
//        cout << "brightness: " << setprecision(2) << linearBrightness << " " << colorBrightnessLinear(color) << " " << targetLinearBrightness << endl;
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

    float squaredEuclidianLinearDistance(Color8 & color1, Color8 & color2) {
        // https://en.wikipedia.org/wiki/Color_difference
        auto c1 = color1.inLinearRGB();
        auto c2 = color2.inLinearRGB();
        
        auto r = (c1[0] + c2[0]) / 2.f;
        auto dR = c1[0] - c2[0];
        auto dG = c1[1] - c2[1];
        auto dB = c1[2] - c2[2];
        
        return dR*dR*(2 + r) + 4*dG*dG + dB*dB*(3 - r);
    }
}
