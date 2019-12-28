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

    }
}
