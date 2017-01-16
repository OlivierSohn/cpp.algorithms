
#include "gtest/gtest.h"

#include "public.h"
// if markov_utils.hpp is included through public.h or private.h, we have a compilation error
// but through he.h, it's fine... weird !!!
#include "he.h"


#include "test_cast.cpp"
#include "test_markov_chain.cpp"
#include "test_cpp.cpp"
#include "test_types.cpp"
#include "test_pool.cpp"
#include "test_sort.cpp"
#include "test_hash_table.cpp"
#include "measure_containers.cpp"
