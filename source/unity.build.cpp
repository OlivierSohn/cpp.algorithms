/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

#include "private.h"

namespace imajuscule {
    namespace fft {

        a64::vector<int8_t> & getFFTTmp() {
            thread_local a64::vector<int8_t> v;
            return v;
        }
    }
}

#include "thread.cpp"
#include "profiling.cpp"
#include "dsp.convolution.cpp"
#include "hash_table.cpp"
#include "pool.adaptive_stack.cpp"
#include "print_type.cpp"
#include "range.cpp"
#include "enum.cpp"
#include "interpolation.cpp"
#include "string.cpp"
#include "sort_utils.cpp"
#include "string.manip.cpp"
#include "debugging.cpp"
#include "scheduler.cpp"
#include "indentedStream.cpp"
#include "color.cpp"
