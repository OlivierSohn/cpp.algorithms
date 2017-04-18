
#include "gtest/gtest.h"

#ifdef MEASURE_PERFS
#error "redefinition of MEASURE_PERFS"
#endif
#define MEASURE_PERFS 0 // to make non regression tests run faster

#include "public.h"

#include <thread>

#include "test_fft_fbins.cpp"
#include "test_fft_signal.cpp"
#include "test_vdsp.cpp"
#include "test_dsp.convolution.cpp"
#include "test_benchmark.cpp"
#include "test_gradient_descent.cpp"
#include "test_peaks.cpp"
#include "test_cast.cpp"
#include "test_cyclic.cpp"
#include "test_fft.cpp"
#include "test_markov_chain.cpp"
#include "test_freelist.cpp"
#include "test_cpp.cpp"
#include "test_types.cpp"
#include "test_allocators.cpp"
#include "test_sort.cpp"
#include "test_hash_table.cpp"

#if MEASURE_PERFS
#  include "measure_containers.cpp"
#endif
