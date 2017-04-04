/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule {
    std::minstd_rand & shuffle_rng_engine() {
        static std::minstd_rand engine;
        return engine;
    }
    
} // NS imajuscule

