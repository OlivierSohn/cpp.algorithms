/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule {
    namespace fft {
        
        template<typename TAG, typename T>
        struct RealInput_;
        
        template<typename TAG, typename T>
        struct Context_;
        
        template<typename TAG, typename T>
        struct ScopedContext_ {
            using CTXT = Context_<TAG, T>;
            ScopedContext_(int lgSize) :
            ctxt(CTXT::create(lgSize))
            {}
            
            ~ScopedContext_() {
                CTXT::destroy(ctxt);
            }
            
            typename CTXT::type ctxt;
            auto get() const { return ctxt; }
        };

    }
}

