
#include "gtest/gtest.h"

#ifdef MEASURE_PERFS
#error "redefinition of MEASURE_PERFS"
#endif
#define MEASURE_PERFS 0 // to make non regression tests run faster

#include "public.h"

#include "test.h"

#include <thread>

#include "test_staticvector_lockfree_scmp.cpp"
#include "test_staticvector_singlethread.cpp"
#include "test_pair_array_distant.cpp"
#include "test_pair_array_local.cpp"
#include "test_fifo.cpp"
#include "test_fifo1.cpp"
#include "test_angles.cpp"
#include "test_locks.cpp"
#include "test_rng.cpp"
#include "test_containers.cpp"
#include "test_math.cpp"
#include "test_fft_fbins.cpp"
#include "test_fft_signal.cpp"
#include "test_vdsp.cpp"
#include "test_dsp.convolution.cpp"
#include "test_dsp.spatialize.cpp"
#include "test_benchmark.cpp"
#include "test_range_search.cpp"
#include "test_global_search.cpp"
#include "test_gradient_descent.cpp"
#include "test_peaks.cpp"
#include "test_cast.cpp"
#include "test_cyclic.cpp"
#include "test_fft.cpp"
#include "test_markov_chain.cpp"
#include "test_scoped.cpp"
#include "test_freelist.cpp"
#include "test_cpp.cpp"
#include "test_types.cpp"
#include "test_allocators.cpp"
#include "test_sort.cpp"
#include "test_hash_table.cpp"

#if MEASURE_PERFS
#  include "measure_containers.cpp"
#  include "test_measure_multiplyadd_contiguity.cpp"
#endif
