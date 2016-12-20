#pragma once

namespace imajuscule {
    
    template<typename T>
    struct Traits;
    
    template<>
    struct Traits<int> {
        static int prehash(int value) {
            return value;
        }
    };

} // NS imajuscule

