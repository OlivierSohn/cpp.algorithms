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
#include "binary.h"
#include "angles.h"
#include "indentedStream.h"
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
#include "profiling.h"
#include "rng.hpp"
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
#include "bounded_lifo.hpp"
#include "hash.hpp"
#include "hash_table.hpp"
#include "merge_sort.hpp"
#include "logger.h"

#include "fft.interface.hpp"

#if __APPLE__
# include "fft.impl.acc.hpp"
#endif

#include "fft.roots.hpp"
#include "fft.impl.imj.hpp"
#include "fft.hpp"

#include "peaks.hpp"
#include "markov_chain.hpp"
#include "dsp.filter.hpp"
#include "dsp.filter.fir.hpp"
#include "dsp.convolution.hpp"
#include "dsp.convolution.finegrained.hpp"
#include "dsp.convolution.combine.hpp"
#include "dsp.spatialize.hpp"
#include "dsp.compresser.hpp"
#include "enum.h"
#include "interpolation.h"
#include "scheduler.h"
#include "scoped.h"
#include "file2string.h"
#include "os.storage.h"
#include "bsonutils.hpp"
#include "utf8.hpp"
#include "bsonparser.hpp"
#include "bsonwriter.hpp"


#if __APPLE__
# include "measure_multiplyadd_contiguity.hpp" // for test on ios
#endif
