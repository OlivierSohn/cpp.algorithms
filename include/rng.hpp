/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule {
    struct rng {
        static std::mt19937 & mersenne() {
            thread_local std::mt19937 mersenne_engine_(std::random_device{}());
            return mersenne_engine_;
        }
    };
    
} // NS imajuscule

