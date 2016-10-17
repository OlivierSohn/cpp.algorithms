#pragma once

namespace imj {

    constexpr bool zero_or_powOf2(int x) {
        return (x & (x-1)) == 0;
    }
    
    static inline bool zero_or_powOf2(size_t x) {
        while(x) {
            if( 1 == x % 2 ) {
                return false;
            }
            x /= 2;
        }
        return true;
    }
    
    
    static size_t pow2(int power) {
        size_t res = 1;
        while(power) {
            res *= 2;
            --power;
        }
        return res;
    }
    
    template<typename T>
    static inline unsigned int log2( T x )
    {
        unsigned int ans = 0 ;

        while( x /= 2 )  {
            ans++;
        }
        return ans ;
    }

} // NS imj

