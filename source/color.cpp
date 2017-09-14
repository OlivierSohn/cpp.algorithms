
namespace imajuscule {
    
    void alphaCompositeColor(Color8 & color, unsigned char overlayed, float overlayAlpha) {
        auto const & c = color.inLinearRGB();
        std::array<unsigned char, 3> linear;
        
        for(int i=0; i<3; ++i) {
            linear[i] = (int)(.5f + (c[i] * (1 - overlayAlpha) + overlayed * overlayAlpha));
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
        rY * c[0] / 255.f +
        gY * c[1] / 255.f +
        bY * c[2] / 255.f;
    }
    
    void adjustColorBrightness(BrightnessAdjustment adj, Color8 & color, float targetLinearBrightness) {
        using namespace std;
        auto linearBrightness = colorBrightnessLinear(color);
        if(!linearBrightness) {
            return;
        }
        
        if(linearBrightness == targetLinearBrightness) {
            return;
        }
        
        if(adj == BrightnessAdjustment::AlphaCompositing) {
            float alpha;
            int overlay;
            if(linearBrightness < targetLinearBrightness) {
                // white is 1
                overlay = 255;
                
                // c is brightness
                // alpha compositing: targetBrightness = (1-alpha) * brightness + alpha * 1
                // alpha * (1 - brightness) = targetBrightness - brightness
                alpha = 1 - ((1 - targetLinearBrightness) / (1-linearBrightness));
            }
            else {
                // black is 0
                overlay = 0;
                
                // c is brightness
                // alpha compositing: targetBrightness = (1-alpha) * brightness
                // alpha * (- brightness) = targetBrightness - brightness
                alpha = 1 - (targetLinearBrightness / linearBrightness);
                
                // note that this branch is equivalent to Scaling method
            }
            assert(alpha <= 1);
            assert(alpha >= 0);
            alphaCompositeColor(color, overlay, alpha);
        }
        else {
            auto ratio = targetLinearBrightness / linearBrightness;
            
            std::array<unsigned char, 3> linearCandidate;
            auto const & c = color.inLinearRGB();
            
            for(int i=0; i<3; ++i) {
                // it is easy to match brightness this way
                // but it might look better to do alpha compositing with white instead
                auto res = (int)(c[i] * ratio);
                if(res > 255) {
                    cout << "color brightness adjustment : saturation" << endl;
                    res = 255;
                }
                linearCandidate[i] = res;
            }
            color.setInLinear(std::move(linearCandidate));
        }
        cout << "brightness: " << setprecision(2) << linearBrightness << " " << colorBrightnessLinear(color) << " " << targetLinearBrightness << endl;
    }
    
    // todo ponderate this with min or max brightness
    // or apply transformation on colors beforehand
    float squaredEuclidianDistance(Color8 & color1, Color8 & color2) {
        // https://en.wikipedia.org/wiki/Color_difference
        auto c1 = color1.inLinearRGB();
        auto c2 = color2.inLinearRGB();
        
        auto r = (c1[0] + c2[0]) / 2.f;
        auto dR = c1[0] - c2[0];
        auto dG = c1[1] - c2[1];
        auto dB = c1[2] - c2[2];
        
        return dR*dR*(2 + r/256) + 4*dG*dG + dB*dB*(2 + (255-r)/256);
    }
}
