/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule
{
    struct Object {
        virtual ~Object() = default;
        Object() = default;
    };

    struct NonCopyable : public Object {
        NonCopyable() =default;
        NonCopyable(const NonCopyable &) = delete;
        NonCopyable & operator=(const NonCopyable&) = delete;
        
        NonCopyable(NonCopyable &&) = default;
        NonCopyable& operator = (NonCopyable &&) = default;
    };
}
