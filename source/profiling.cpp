/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule
{
    namespace profiling {
        void pollute_cache()
        {
            std::vector<unsigned int> pollution(100000);
            unsigned int i = 0;
            for(auto & p : pollution) {
                p = i++;
            }
            i = 0;
            for(auto & p : pollution) {
                i += p;
            }
            std::cout << ((i%2)?"." : ":");
        };        
    }
}
