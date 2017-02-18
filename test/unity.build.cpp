
#include "gtest/gtest.h"

#ifdef MEASURE_PERFS
#error "redefinition of MEASURE_PERFS"
#endif
#define MEASURE_PERFS 0 // to make non regression tests run faster

#include "public.h"

#include <thread>

// if markov_utils.hpp is included through public.h or private.h, we have a compilation error
// but through he.h, it's fine... weird !!!


#include "test_cast.cpp"
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
