/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

#pragma once

#include <algorithm>
#include <array>
#include <atomic>
#include <cassert>
#include <cctype>
#include <chrono>
#include <cstddef>
#ifdef IMJ_LOG_MEMORY
#  include <cstdlib>
#endif
#include <cstring>
#include <ctime>
#include <cxxabi.h>
#include <fstream>
#include <functional>
#include <initializer_list>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <list>
#include <locale>
#include <map>
#include <memory>
#include <mutex>
#include <new>
#include <numeric>
#include <random>
#include <set>
#include <sstream>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <string>
#include <thread>
#include <type_traits>
#include <typeinfo>
#include <unordered_map>
#include <utility>
#include <variant>
#include <vector>

#ifdef __APPLE__
#include "TargetConditionals.h"
#endif

#if defined (__unix__) || (defined (__APPLE__) && defined (__MACH__))
# include <unistd.h>
#endif

#if (defined (__APPLE__) && defined (__MACH__))

# include <mach/thread_policy.h>

# if TARGET_OS_IPHONE
# else
#  include <sys/sysctl.h>
# endif

#endif

#if IMJ_WITH_OPENCL
#  ifdef __APPLE__
#    include <OpenCL/opencl.h>
#  else
#    include <CL/cl.h>
#  endif // __APPLE__
#endif // WITH_OPENCL

#if __APPLE__
# include <Accelerate/Accelerate.h>
# include "api.accelerate.h"
#endif

#if __has_include(<optional>)
#   include <optional>
#elif __has_include(<experimental/optional>)
#   include <experimental/optional>
#else
#   error Must have an optional type, either from <optional> or if not supported from <experimental/optional>.
#endif

#include "thirdparty/atomic_queue/atomic_queue.h"

#include "array.utils.hpp"
#include "csv.writer.hpp"
#include "log.stack.h"
#include "maybe.atomic.hpp"
#include "likely.h"
#include "spinlock.h"
#include "log.h"
#include "imj.assert.h"
#include "split.hpp"
#include "cconstarray.hpp"
#include "convex_hull.h"
#include "color.h"
#include "matrix.hpp"
#include "linalg.h"
#include "binary.h"
#include "angles.h"
#include "indentedStream.h"
#include "either.hpp"
#include "logging.h"
#include "container.hpp"
#include "optional.h"
#include "fifo1.hpp"
#include "object.h"
#include "safe_cast.hpp"
#include "gen.names.h"
#include "available_indexes.hpp"
#include "fifo.hpp"
#include "slidingsum.h"
#include "meta.hpp"
#include "numtraits.h"
#include "math.hpp"
#include "flt_math.hpp"
#include "freelist.hpp"
#include "complex.hpp"
#include "range.h"
#include "strplot.h"
#include "thread.h"
#include "rng.hpp"
#include "profiling.h"
#include "optimization.hpp"
#include "optimization.global.hpp"
#include "optimization.local.hpp"
#include "optimization.local.ranged.hpp"
#include "debugging.h"
#include "print_type.hpp"
#include "pool.adaptive_stack.h"
#include "allocator.adaptive_stack.hpp"
#include "allocator.aligned.hpp"
#include "allocated_containers.h"
#include "cyclic.h"
#include "slidingaverage.h"
#include "stackvector.hpp"
#include "pair.array.hpp"
#include "staticvector.lockfree.scmp.hpp"
#include "staticvector.singlethread.hpp"
#include "staticvector.hpp"
#include "fifo.lockfree.scmp.hpp"
#include "forward.list.lockfree.scmp.hpp"
#include "string.manip.h"
#include "iter_range.hpp"
#include "sort_utils.hpp"
#include "insertion_sort.hpp"
#include "heap_sort.hpp"
#include "n_insertion_sort.hpp"
#include "bounded_lifo.hpp"
#include "hash.hpp"
#include "hash_table.hpp"
#include "merge_sort.hpp"
#include "logger.h"

#include "math.roots.hpp"

#include "fft.interface.hpp"
#include "fft.costs.hpp"

#if __APPLE__
# include "fft.impl.acc.hpp"
#endif

#include "fft.roots.hpp"
#include "fft.impl.imj.hpp"
#include "fft.hpp"

#include "smooth.hpp"
#include "peaks.hpp"
#include "markov_chain.hpp"

extern "C"
{
#include "c.h"
}

#include "os.storage.h"
#include "samples.h"
#include "dsp.resample.h"
#include "read.wav.h"

#include "dsp.epsilon.hpp"
#include "dsp.cost.hpp"
#include "dsp.host.simulation.hpp"
#include "dsp.convolution.optimization.hpp"
#include "dsp.filter.hpp"
#include "dsp.filter.fir.hpp"
#include "dsp.subsampled.hpp"
#include "dsp.convolution.async.hpp"
#include "dsp.convolution.hpp"
#include "dsp.convolution.finegrained.hpp"
#include "dsp.delayed.hpp"
#include "dsp.convolution.scale.hpp"
#include "dsp.convolution.scale.custom.hpp"
#include "dsp.convolution.scaling.hpp"
#include "dsp.convolution.split.hpp"
#include "dsp.convolution.combine.hpp"

#if IMJ_WITH_OPENCL
#  include "dsp.convolution.gpu.hpp"
#endif // WITH_OPENCL

#include "dsp2.hpp"
#include "dsp2.filter.fir.hpp"
#include "dsp2.convolution.hpp"
#include "dsp2.convolution.split.hpp"
#include "dsp2.convolution.scale.custom.hpp"

#include "dsp.impulseresponse.hpp"
#include "dsp.spatialize.hpp"
#include "dsp.reverb.hpp"
#include "dsp.reverbs.hpp"
#include "dsp.convolution.bycbsize.hpp"
#include "dsp.compresser.hpp"
#include "enum.h"
#include "interpolation.h"
#include "scheduler.h"
#include "scoped.h"
#include "file2string.h"
#include "bsonutils.hpp"
#include "utf8.hpp"
#include "bsonparser.hpp"
#include "bsonwriter.hpp"


#if __APPLE__
# include "measure_multiplyadd_contiguity.hpp" // for test on ios
#endif
