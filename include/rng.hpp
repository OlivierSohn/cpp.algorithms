/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule {
    
    enum class SEEDED {
        Yes,
        No
    };
    
    template<SEEDED S>
    std::mt19937 & mersenne() {
        thread_local std::mt19937 mersenne_engine_(std::random_device{}());
        return mersenne_engine_;
    }


template<SEEDED S>
std::ranlux24_base & lagged_fibonacci() {
    thread_local std::ranlux24_base e;
    return e;
}
    
} // NS imajuscule

