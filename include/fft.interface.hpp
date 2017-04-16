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
        struct RealOutput_;
        
        template<typename TAG, typename T>
        struct Context_;
        
        template<typename TAG, typename T>
        struct ScopedContext_ {
            using CTXT = Context_<TAG, T>;
            ScopedContext_(int size) :
            ctxt(CTXT::create(size) )
            {}
            
            ~ScopedContext_() {
                CTXT::destroy(ctxt);
            }
            
            typename CTXT::type ctxt;
            auto get() const { return ctxt; }
        };
        
        template<typename TAG, typename T>
        struct Algo_;
        
        namespace slow_debug {
            template<typename TAG, typename CONTAINER>
            struct UnwrapFrequencies;
            
            template<typename TAG, typename CONTAINER>
            struct UnwrapSignal;
            
            template<typename TAG, typename CONTAINER>
            auto unwrap_frequencies(CONTAINER const & c, int size) {
                UnwrapFrequencies<TAG, CONTAINER> u;
                return u.run(c, size);
            }
            
            template<typename TAG, typename CONTAINER>
            auto unwrap_signal(CONTAINER const & c, int size) {
                UnwrapSignal<TAG, CONTAINER> u;
                return u.run(c, size);
            }
        } // NS slow_debug
    } // NS fft
} // NS imajuscule

