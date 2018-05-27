/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule
{
    struct Logger {
        template <class... Args>
        static void err(Args&&... args) {
            LG(ERR, string_format(std::forward<Args>(args)...).c_str());
        }
    };

}
