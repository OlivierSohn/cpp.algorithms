/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule
{
    namespace scoped {
        
        // todo support multithreading
        template<typename Parent>
        struct OnLeavingLast : public Parent {
            using Parent::f;
            using Parent::n;
            
            OnLeavingLast() { ++n; }
            
            ~OnLeavingLast() {
                auto new_value = --n;
                if(0 == new_value) {
                    // there is no other instance in a parent scope
                    f();
                }
                else {
                    // an instance exists in a parent scope
                }
            }
        };

    }
}
