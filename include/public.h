/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

#pragma once

#include <array>
#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstring>
#include <cstddef>
#include <ctime>
#include <cxxabi.h>
#include <functional>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <list>
#include <map>
#include <new>
#include <unordered_map>
#include <memory>
#include <random>
#include <set>
#include <sstream>
#include <stdio.h>
#include <string>
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

#include <mach/thread_policy.h>

# if TARGET_OS_IPHONE
// those arre commented out in the header on ios
kern_return_t thread_policy_set
(
	thread_act_t thread,
	thread_policy_flavor_t flavor,
	thread_policy_t policy_info,
	mach_msg_type_number_t policy_infoCnt
 );
kern_return_t thread_policy_get
(
	thread_act_t thread,
	thread_policy_flavor_t flavor,
	thread_policy_t policy_info,
	mach_msg_type_number_t *policy_infoCnt,
	boolean_t *get_default
 );
# else
#  include <sys/sysctl.h>
# endif

#endif

#if __APPLE__
# include <Accelerate/Accelerate.h>
#endif

#if __APPLE__
# include "api.accelerate.h"
#endif

#include "defines.h"
#include "logging.h"
#include "object.h"
#include "safe_cast.hpp"
#include "available_indexes.hpp"
#include "cyclic.h"
#include "slidingsum.h"
#include "slidingaverage.h"
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
#include "static_vector.h"
#include "string.manip.h"
#include "iter_range.hpp"
#include "sort_utils.hpp"
#include "insertion_sort.hpp"
#include "heap_sort.hpp"
#include "bounded_lifo.hpp"
#include "hash.hpp"
#include "hash_table.hpp"
#include "merge_sort.hpp"

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
#include "dsp.convolution.hpp"
#include "dsp.convolution.finegrained.hpp"
#include "dsp.convolution.benchmarks.hpp"
#include "enum.h"
#include "interpolation.h"
#include "scheduler.h"
