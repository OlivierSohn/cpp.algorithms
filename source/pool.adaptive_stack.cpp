
namespace imajuscule {
#ifndef NDEBUG
    AdaptiveStack::State thread_local AdaptiveStack::state = AdaptiveStack::Growing;
#endif

    void AdaptiveStack::allocate_overflow() {
        overflow.reset( new AdaptiveStack(buffer.size()) );
    }
    
    AdaptiveStack & AdaptiveStack::getThreadLocalInstance() {
        thread_local AdaptiveStack instance;
        return instance;
    }
    
} // ns imajuscule
