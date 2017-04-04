/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule
{
    struct tm * getTime() {
        time_t result;
        result = time(nullptr);
        
        struct tm * pTime = nullptr;
#ifdef _WIN32
        static struct tm time;
        pTime = &time;
        localtime_s(pTime, &result);
#else
        pTime = localtime(&result);
#endif
        return pTime;
    }
}
