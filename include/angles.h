/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule {
    static inline float modulo_angle(float a) {
        auto res = a + M_PI;
        
        auto r = res / (2*M_PI);
        if(r < 0) {
            res += (2*M_PI) * (1 + (int)-r);
        }
        else {
            res -= (2*M_PI) * (int)r;
        }
        
        return res - M_PI;
    }
}
