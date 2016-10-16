/*
* Tests of imj::merge_sort
*/

#include <stdio.h>

#include <algorithm>
#include <vector>
#include <list>
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
                performance_test(length);
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
            
            std::srand(0);

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
        void performance_test(int size, int n=3, bool compare = true) {
            auto imj_time = 0.f;
            auto std_time = 0.f;

            auto c = make_container_no_repeat(size);

            for(auto i=0; i<n; i++) {
                Shuffle(c);
                
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
            
            if(compare) {
                printf("\nimj: %.4fs\n", imj_time);
                printf("\nstd: %.4fs\n", std_time);                
            }
        }
        
    };
    
    template< typename Container >
    void testMergeSort() {
        
        cout << "---" << endl;
        cout << "test MergeSort for" << endl;
        PRINT_TYPE(Container);
        cout << "---" << endl;
        
        std::vector<AlgoType> algo_types{
            RECURSIVE,
            SEQUENTIAL,
            SEQUENTIAL_CACHE_OPTIMIZED
        };
        
        for(auto t : algo_types) {
            Test<Container> test;
            
            test.length_power = 4;
            test.length_power_performance = 4;
            
            std::cout << "mode " << ((t==RECURSIVE)? "RECURSIVE" : ((t==SEQUENTIAL)?"SEQUENTIAL":"SEQUENTIAL_CACHE_OPTIMIZED")) << std::endl;
            
            test.setSorter([t](Container & c) {
                merge_sort(c.begin(), c.end(), t);
            });
            
            test.run();
        }
    }
    
    template< typename Container >
    void testInsertionSort() {
        cout << "---" << endl;
        cout << "test InsertionSort for" << endl;
        PRINT_TYPE(Container);
        cout << "---" << endl;
        
        Test<Container> test;
        test.length_power = 2;
        test.length_power_performance = 2;
        
        test.setSorter([](Container & c) {
            insertion_sort(c.begin(), c.end());
        });
        
        test.run();
    }
    
    
} // NS imj

using namespace imj;

TEST(Algorithm, MergeSort) {
    // disable printf buffering
    setbuf(stdout, NULL);
    
    testMergeSort< vector<int> >();
    testMergeSort< list<int> >();
}

/*
TEST(Algorithm, MergeSort_profile) {
    // disable printf buffering
    setbuf(stdout, NULL);
    
    imj::Test<vector<int>> test;
    test.setAlgoType(SEQUENTIAL);
    test.performance_test(1000000, 100, false);
}
*/


TEST(Algorithm, InsertionSort) {
    // disable printf buffering
    setbuf(stdout, NULL);
    
    testInsertionSort< vector<int> >();
    testInsertionSort< list<int> >();
}
