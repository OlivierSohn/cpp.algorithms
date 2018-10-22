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
#include "dsp.subsampled.hpp"
#include "dsp.convolution.hpp"
#include "dsp.convolution.finegrained.hpp"
#include "dsp.delayed.hpp"
#include "dsp.convolution.combine.hpp"


namespace imajuscule
{
  template<typename C>
  void applySetup(C&c, typename C::SetupParam const & p) {
    c.applySetup(p);
  }
  
  template<typename C>
  void applySetup(Delayed<C>&c, typename Delayed<C>::SetupParam const & p) {
    c.setTheDelay(p.delay);
    applySetup(c.getInner(),p.innerParams);
  }
  template<LatencySemantic L, typename C>
  void applySetup(SubSampled<L, C>&c, typename C::SetupParam const & p) {
    applySetup(c.getInner(), p);
  }
  
  template<typename A, typename B>
  void applySetup(SplitConvolution<A,B> &c, typename SplitConvolution<A,B>::SetupParam const & p) {
    applySetup(c.getA(), p.aParams);
    applySetup(c.getB(), p.bParams);
    if(!c.getB().isValid()) { // handle case where lower resolution tail is not used
      c.setSplit(noSplit);
    }
    else {
      c.setSplit( c.getB().getLatency() - c.getA().getLatency() );
    }
  }
  template<typename A, typename B>
  void applySetup(SplitConvolution<A,ScaleConvolution<B>> &c,
                  typename SplitConvolution<A,ScaleConvolution<B>>::SetupParam const & p) {
    applySetup(c.getA(), p.aParams);
    applySetup(c.getB(), p.bParams);
    c.setSplit(pow2(c.getB().getEarlyDroppedConvolutions()) - 1);
  }
  
  template<typename SP, typename T, typename FFTTag>
  void prepare(SP const & params,
               FinegrainedPartitionnedFFTConvolution<T,FFTTag> & rev,
               int const n_scales,
               int const scale_sz ) {    
    assert(n_scales == 1);
    applySetup(rev, params);
  }
  
  template<typename SP, typename T, typename FFTTag>
  void prepare(SP const & params,
               ZeroLatencyScaledFineGrainedPartitionnedConvolution<T,FFTTag> & rev,
               int const n_scales,
               int const scale_sz ) {
    // split = latB - latA
    // scale_sz = (1 + 2*(delay + latA)) - latA
    // scale_sz - 1 = 2*delay + latA
    // delay = 0.5 * (scale_sz - 1 - latA)
    // delay = 0.5 * (scale_sz - 1 - (2*partition_sz-1))
    // delay = (0.5 * scale_sz) - partition_sz
    assert(0 == scale_sz % 2);
    int const delay = (scale_sz / 2) - params.partition_size;
    assert(delay >= 0);
    
    // set the delays
    
    // we disable the unused scales by setting the partition size to 0.
    auto zero = SP::makeInactive();
    std::array<SP, 4> ps {
      params,
      (n_scales >= 2)?params:zero,
      (n_scales >= 3)?params:zero,
      (n_scales >= 4)?params:zero
    };
    assert(n_scales <= 4);
    applySetup(rev,
               {
                 {},
                 {
                   ps[0],
                   {
                     (n_scales >= 2)?delay:0,
                     {
                       ps[1],
                       {
                         (n_scales >= 3)?delay:0,
                         {
                           ps[2],
                           {
                             (n_scales >= 4)?delay:0,
                             ps[3]
                           }
                         }
                       }
                     }
                   }
                 }
               }
               );
  }
}

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
