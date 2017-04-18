/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule {
    namespace fft {
        
        template<typename TAG, typename T>
        struct RealSignal_;
        
        template<typename TAG, typename T>
        struct RealFBins_;
        
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
        struct Contexts_ {
            using Context  = Context_<TAG, T>;
            using ContextT = typename Context::type;
            
            static Contexts_ & getInstance() {
                // ok to have static variable in header because class is templated
                // (cf. test ThreadLocal)
                thread_local Contexts_ ctxt;
                
                return ctxt;
            }
            
            ContextT getBySize(int size) {
                assert(size > 0);
                assert(is_power_of_two(size));
                auto index = power_of_two_exponent(size);
                if(index >= contexts.size()) {
                    contexts.resize(index+1);
                }
                auto & ret = contexts[index];
                if(!ret) {
                    ret = Context::create(size);
                }
                return ret;
            }
        private:
            Contexts_() {
                contexts.resize(20);
            }
            std::vector<ContextT> contexts;
        };
        
        template<typename TAG, typename T>
        struct Algo_;
        
        namespace slow_debug {
            template<typename TAG, typename CONTAINER>
            struct UnwrapFrequenciesRealFBins;
            
            template<typename TAG, typename CONTAINER>
            struct UnwrapSignal;
            
            template<typename TAG, typename CONTAINER>
            auto unwrap_frequencies(CONTAINER const & c, int size) {
                UnwrapFrequenciesRealFBins<TAG, CONTAINER> u;
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

