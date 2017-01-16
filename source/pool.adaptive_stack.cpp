
namespace imajuscule {
    
    AdaptiveStack * AdaptiveStack::instance = nullptr;
#ifndef NDEBUG
    AdaptiveStack::State AdaptiveStack::state = AdaptiveStack::Growing;
#endif
    
    AdaptiveStack & AdaptiveStack::getInstance() {
        if(!instance) {
            instance = new AdaptiveStack();
        }
        return *instance;
    }
    
    void AdaptiveStack::allocate_overflow() {
        overflow.reset( new AdaptiveStack(buffer.size()) );
    }
    
} // ns imajuscule
