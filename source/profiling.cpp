/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule
{
    namespace profiling {
      char progressChar(unsigned int i) {
        switch(i%4) {
          case 0:
          return '-';
          case 1:
          return '\\';
          case 2:
          return '|';
          case 3:
          return '/';
          default:
          assert(0);
        }
        return ' ';
      }

      std::string progressStr(int sz, unsigned int i) {
        std::string s(sz, '.');
        s[0] = '[';
        s[sz-1] = ']';
        s[1+i%(sz-2)] = 'O';
        return s;
      }

        void pollute_cache()
        {
          // to have alternating characters we need some state:
          static unsigned int n = 0;
          ++n;

            std::vector<unsigned int> pollution(100000);
            unsigned int i = 0;
            for(auto & p : pollution) {
                p = i++;
            }
            i = 0;
            for(auto & p : pollution) {
                i += p;
            }
            // erase last char after the flush so that it is visible until a new flush occurs.
            constexpr int szProgress = 20;
            std::cout
            << progressStr(szProgress, i+n)
            << std::flush
            << std::string(szProgress, '\b')
            << std::string(szProgress, ' ')
            << std::string(szProgress, '\b');
        };
    }
}
