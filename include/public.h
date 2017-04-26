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

#if __APPLE__
# include <Accelerate/Accelerate.h>
#endif

#if __APPLE__
# include "api.accelerate.h"
#endif

#include "profiling.h"
#include "defines.h"
#include "logging.h"
#include "object.h"
#include "safe_cast.hpp"
#include "available_indexes.hpp"
#include "cyclic.h"
#include "slidingaverage.h"
#include "meta.hpp"
#include "numtraits.h"
#include "math.hpp"
#include "flt_math.hpp"
#include "freelist.hpp"
#include "complex.hpp"
#include "range.h"
#include "strplot.h"
#include "rng.hpp"
#include "optimization.hpp"
#include "optimization.global.hpp"
#include "optimization.local.hpp"
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
