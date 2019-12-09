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
namespace audio {
std::ostream & operator << (std::ostream & os, const ResampleSincStats& s) {
    os << "resample    : ";
    if (s.dt_resample) {
        os << s.dt_resample->count()/1000000.f << " s";
    }
    else {
        os << "?";
    }
    os << std::endl;
    os << "read source : ";
    if (s.dt_read_source) {
        os << s.dt_read_source->count()/1000000.f << " s";
    }
    else {
        os << "?";
    }
    os << std::endl;
    return os;
}
}
}

#include "memory.cpp"
#include "log.cpp"
#include "log.stack.cpp"
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
#include "IndentedStream.cpp"
#include "color.cpp"
#include "gen.names.cpp"
#include "os.storage.cpp"
#include "file2string.cpp"
#include "read.wav.cpp"
#include "write.wav.cpp"
#include "samples.cpp"
