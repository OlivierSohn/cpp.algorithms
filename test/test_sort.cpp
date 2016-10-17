/*
* Tests of imj::merge_sort
*/

#include <stdio.h>

#include <algorithm>
#include <vector>
#include <list>
#include <functional>
#include <ctime>

#include "gtest/gtest.h"

#include "sort_utils.hpp"
#include "print_type.hpp"

#include "merge_sort.hpp"
#include "insertion_sort.hpp"

using namespace std;
using namespace imj;

namespace imj {
    
    template< typename Container >
    struct Test {

        int length_power = 0;
        int length_power_performance = 0;

        Test() {
            // flush logs at each print
            setbuf(stdout, NULL);
        }
        
        void run() {
            run_performance();
            run_logic();
        }
        
        void run_logic()
        {
            test_is_sorted();
                
            test_small_containers();
            
            for(auto size = 0; size <= 4; ++size) {
                test_permutations(size);
            }

            auto length = 1;
            for(int i=0; i<=length_power; i++) {
                test_n_permutations( length, 10 );
                length *= 10;
            }
        }
    
        void run_performance() {
            auto length = 1;
            for(int i=0; i<=length_power_performance; i++) {
                if(i < length_power_performance - 1) {
                    length *= 10;
                    continue;
                }
                Container ref;
                for(int j=0; j<2; j++) {
                    ref = performance_test(length, 3, true, j==0?0:&ref);
                }
                length *= 10;
            }
        }
        
        void setSorter(std::function<void(Container&)> f) {
            sorter = f;
        }
    private:
        std::function<void(Container &)> sorter;
        
        using Containers = vector<Container>;
                
        auto make_container_with_repeats(int size, int n_repeats, int index_repeat) 
        {
            auto c = make_container_no_repeat(size);
            
            auto it = c.begin();
            advance(it, index_repeat);
            auto repeat_value = *it;
            ++it;
            
            int count = 0;
            for(; it != c.end(); ++it) {
                if(++count >= n_repeats) {
                    break;
                }
                *it = repeat_value;
            }
            
            it = c.begin();
            for(; it != c.end(); ++it) {
                if(++count >= n_repeats) {
                    break;
                }
                *it = repeat_value;
            }
            return c;
        }
        
        auto make_container_no_repeat(int size) 
        {
            Container c(size);
            iota(c.begin(), c.end(), 0);
            return c;
        }
        
        auto make_container_many_repeats(int size, int number_different_values = 3) 
        {
            Container c(size);
            for(auto & v : c) {
                v = std::rand() % number_different_values;
            }
            return c;
        }
        
        void check_sort(Container const & unsorted_const) {
            auto unsorted = unsorted_const;
            bool was_sorted = IsSorted(unsorted);
            Container unsorted_copy = unsorted;
            ASSERT_EQ(unsorted_copy, unsorted);
           
            sorter(unsorted);
            if(!was_sorted) {
                ASSERT_NE( unsorted_copy, unsorted);
            }
            
            StdSort(unsorted_copy);
            ASSERT_EQ(unsorted_copy, unsorted) << logstr(unsorted_const);
        }

        void test_permutations(Container && v) {

            // The container needs to be sorted at the first call of next_permutation
            // so that subsequent calls traverse all possible permutations.

            StdSort(v);
            
            do {                
                check_sort(v);
            }
            while (next_permutation(v.begin(), v.end()));
        }

        void test_permutations( int size )
        {
            printf("sort containers of size %d [all permutations]", size);
            {
                auto c = make_container_no_repeat(size);
                test_permutations(move(c));
            }
                    
            printf(" [with one repeated value] ");
                    
            for(int index_repeat = 0; index_repeat < size; index_repeat++) { 
                for(int n_repeats = 2; n_repeats <= size; ++n_repeats) {
                    auto c = make_container_with_repeats(size, n_repeats, index_repeat);        
                    test_permutations(move(c));
                }
            }

            printf(" [with several repeated values] ");
            std::srand(0);
            for(int i=0; i<100; i++) {
                auto c = make_container_many_repeats(size, 2);
                test_permutations(move(c));
            }
            
            printf("\n");
        }
            
        void test_n_permutations(int size, int n_permutations = 10)
        {
            auto v = make_container_no_repeat(size);
            printf("sort containers of size %d [%d permutation(s)]", size, n_permutations);
            
            ShuffleSeed(0);

            for(int i=0; i<n_permutations; ++i) {
                Shuffle(v);
                check_sort(v);
            }
            printf("\n");
        }

        void test_small_containers() {
            Containers containers{
                {},
                {0},
                {0,0},
                {0,1},
                {1,0},
                {0,1,2},
                {2,1,0},
                {0,2,1},
                {1,2,0},
                {1,2,3,0,4,5},
            };
            
            for(auto & c : containers) {
                check_sort(c);    
            }
        }
         
        void test_is_sorted() {
            Containers sorted {
                {},
                {0},
                {0,0},
                {0,1},
                {0,0,1,2,2,3}
            };
            
            Containers not_sorted {
                {1,0},
                {1,1,0},
                {0,1,0},
                {0,2,1},
                {0,0,2,2,1},
                {0,0,2,2,1,3}
            };
            
            for(auto & c : sorted) {
                ASSERT_TRUE( IsSorted(c) );    
            }
            for(auto & c : not_sorted) {
                ASSERT_FALSE( IsSorted(c) );    
            }
        }        
        
    public:
        auto performance_test(int size, int n=3, bool compare = true, Container * ref = 0) {
            auto imj_time = 0.f;
            auto std_time = 0.f;

            auto c = make_container_no_repeat(size);
            
            ShuffleSeed(0);
            for(auto i=0; i<n; i++) {
                Shuffle(c);
                
                // to reduce measure error
                for(auto repeats = 0; repeats <10; repeats ++)
                {
                    {
                        auto copy = c;
                        clock_t t = clock();
                        sorter(copy);
                        imj_time += (float)(clock() - t)/CLOCKS_PER_SEC;
                    }
                    
                    if(compare)
                    {
                        auto copy = c;
                        clock_t t = clock();
                        StdSort(copy);
                        std_time += (float)(clock() - t)/CLOCKS_PER_SEC;
                    }
                }
            }
            
            if(compare) {
                printf("  %d:\t %.4f (%.4f / %.4f)\n", size, imj_time/std_time, imj_time, std_time);
            }
            
            // to test that seed actually works as intended
            if(ref) {
                if(*ref != c) {
                    std::cout << "!!! ref doesn't match" << std::endl;
                }
            }
            
            return c;
        }
        
    };
    
    template< typename Container, int insertion_sort_below_size >
    void testMergeSort(std::vector<AlgoType> const & algo_types, bool perf_only = false) {
        using iterator = typename Container::iterator;
        
        cout << "---" << endl;
        cout << "MergeSort (" << insertion_sort_below_size << ")" << endl;
        PRINT_TYPE(Container);
        cout << "---" << endl;
        
        for(auto t : algo_types) {
            Test<Container> test;
            
            test.length_power = 2;
            test.length_power_performance = 6;
            
            std::cout << ((t==RECURSIVE)? "RECURSIVE" : ((t==SEQUENTIAL)?"SEQUENTIAL":"SEQUENTIAL_CACHE_OPTIMIZED")) << std::endl;
            
            test.setSorter([t](Container & c) {
                merge_sort<iterator, insertion_sort_below_size>(c.begin(), c.end(), t);
            });
            
            if(perf_only) {
                test.run_performance();
            }
            else {
                test.run_logic();
            }
        }
        
    }
    
    template< typename Container >
    void testInsertionSort() {
        cout << "---" << endl;
        cout << "InsertionSort" << endl;
        PRINT_TYPE(Container);
        cout << "---" << endl;
        
        Test<Container> test;
        test.length_power = 2;
        
        test.setSorter([](Container & c) {
            insertion_sort( c.begin(), c.end() );
        });
        
        test.run_logic();
    }
    
    
} // NS imj

using namespace imj;

#define PERFS 0

#if 0 == PERFS
TEST(Algorithm, MergeSort) {
    std::vector<AlgoType> algo_types{
        RECURSIVE,
        SEQUENTIAL,
        SEQUENTIAL_CACHE_OPTIMIZED
    };
    
    testMergeSort< vector<int>, 0 >(algo_types);
    testMergeSort< vector<int>, 2 >(algo_types);
    testMergeSort< vector<int>, 4 >(algo_types);
    testMergeSort< list<int>, 0>(algo_types);
    testMergeSort< list<int>, 2>(algo_types);
    testMergeSort< list<int>, 4>(algo_types);
}

TEST(Algorithm, InsertionSort) {
    testInsertionSort< vector<int> >();
    testInsertionSort< list<int> >();
}

#else

TEST(Algorithm, MergeSort_profile) {
    std::vector<AlgoType> algo_types{
        SEQUENTIAL
    };

    testMergeSort< vector<int>, 0 >(algo_types, true);
    testMergeSort< vector<int>, 2 >(algo_types, true);
    testMergeSort< vector<int>, 4 >(algo_types, true);
    testMergeSort< vector<int>, 8 >(algo_types, true);
    testMergeSort< vector<int>, 16 >(algo_types, true);
    testMergeSort< vector<int>, 32 >(algo_types, true); // best (for my data / laptop architectures)
    testMergeSort< vector<int>, 64 >(algo_types, true);
    testMergeSort< vector<int>, 128 >(algo_types, true);
    testMergeSort< vector<int>, 256 >(algo_types, true);
}

#endif // PERFS


