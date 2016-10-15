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


using namespace std;
using namespace imj;

namespace imj {
    
    template< typename Container >
    struct Test {

        void run_logic()
        {
            test_is_sorted();
                
            test_small_containers();
            
            for(auto size = 0; size <= 4; ++size) {
                test_permutations(size);
            }

            for(auto multiple = 1; multiple < 100; multiple *= 3) {
                auto size = 100 * multiple;
                
                test_n_permutations( size, 10 );
            }
        }
    
        void run_performance() {
            performance_test(100);
            performance_test(10000);
            performance_test(1000000);
        }
    
    private:
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
           
            imj::merge_sort(unsorted.begin(), unsorted.end());
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
            
        void performance_test(int size) {
            auto imj_time = 0.f;
            auto std_time = 0.f;

            for(auto i=0; i<3; i++) {
                auto c = make_container_no_repeat(size);
                Shuffle(c);
                
                {
                    auto copy = c;
                    clock_t t = clock();                    
                    merge_sort(copy.begin(), copy.end());
                    imj_time += (float)(clock() - t)/CLOCKS_PER_SEC;
                }

                {
                    auto copy = c;
                    clock_t t = clock();
                    StdSort(copy);
                    std_time += (float)(clock() - t)/CLOCKS_PER_SEC;          
                }
            }            
            
            printf("\nimj: %.4fs\n", imj_time);
            printf("\nstd: %.4fs\n", std_time);
        }
        
    };
    
    template< typename Container >
    void test() {
        cout << "---" << endl;
        PRINT_TYPE(Container);
        cout << "---" << endl;
        
        Test<Container> test;
        test.run_performance();
        test.run_logic();
    }
    
    
} // NS imj

TEST(Algorithm, MergeSort) {
    // disable printf buffering
    setbuf(stdout, NULL);
        
    test< vector<int> >();
    //test< vector<float> >();
    test< list<int> >();
    //test< list<float> >();
}

