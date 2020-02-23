
namespace imajuscule {
    namespace fft {
        
        template<typename TAG, typename T>
        struct RealSignal_;
        
        template<typename TAG, typename T, template<typename> typename Allocator>
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
                static Contexts_ ctxt;
                
                return ctxt;
            }
            
            ContextT getBySize(int size) {
                std::lock_guard l(mut);
                
                assert(size > 0);
                assert(is_power_of_two(size));
                auto index = power_of_two_exponent(size);
                if(index >= contexts.size()) {
                    contexts.resize(2*index);
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
          ~Contexts_() {
            for(auto const &c:contexts) {
              if(c) {
                Context::destroy(c);
              }
            }
          }
            std::vector<ContextT> contexts;
            std::mutex mut;
          
          Contexts_(const Contexts_&) = delete;
          Contexts_(Contexts_&&) = delete;
          Contexts_& operator=(const Contexts_&) = delete;
          Contexts_& operator=(Contexts_&&) = delete;
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


/*
 For an intuitive presentation of the 2 methods, see:
 https://www.slideshare.net/GourabGhosh4/overlap-add-overlap-savedigital-signal-processing
 
 For our use-case, overlap save has several advantages:
 . it allows to perform ffts of x _concurrently_ using a single x buffer because we don't need 0-padding
 . it makes less (half) memory writes to y (because only the second half of the spectrum is used)
 . it uses less (half) memory for y (for the same reason mentionned above)
 
 The only inconvenient overlap save has is that:
 . with overlap add, we can do the x 0-padding incrementally, so the worst (peak) number of reads / writes to x is
   maxFftSize/4
 . with overlap save, copying x needs to happen at once, so the worst (peak) number of reads / writes to x is
   maxFftSize/2
 
 When the biggest fft of x and the write of the results to y happen in the same step,
 this inconvenient is entirely compensated by the fact that
 there are less memory writes to y.
 
 But when the biggest fft of x and the write to y happen in _different_ steps
 (like in the finegrained case), this inconvenient leads to a higher worst case step time.
*/
enum class Overlap {
    Add,
    Save
};

constexpr auto overlapMode = Overlap::Save;

} // NS imajuscule

