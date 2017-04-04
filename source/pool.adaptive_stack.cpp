/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

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
