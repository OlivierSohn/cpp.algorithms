/*
* Tests of HashTable
*/

#include <stdio.h>

#include <algorithm>
#include <vector>
#include <list>
#include <functional>
#include <ctime>

#include "gtest/gtest.h"

#include "hash_table.hpp"
#include "print_type.hpp"

using namespace std;
using namespace imj;

namespace imj {
namespace test {
namespace hastable {
    
    template< typename Value >
    struct Test {
        
        void run_logic()
        {
            HashTable<Value> h;
            
            EXPECT_TRUE(h.isConsistent());
            
            EXPECT_TRUE( h.empty() );
            EXPECT_EQ( 0, h.size());
            EXPECT_EQ( 0, h.max_chain_length());
            
            int max_ = 1000;
            
            for(int i=0; i < max_; i++) {
                ASSERT_FALSE(h.has(i));
                h.insert(i);
                ASSERT_TRUE(h.has(i));
                EXPECT_TRUE(h.isConsistent());
                EXPECT_LE(1, h.max_chain_length());
            }
            
            EXPECT_FALSE( h.empty() );
            EXPECT_EQ( max_, h.size());

            for(int i=0; i< max_/2; i++) {
                ASSERT_TRUE(h.has(i));
                h.remove(i);
                ASSERT_FALSE(h.has(i));
                EXPECT_TRUE(h.isConsistent());
                EXPECT_LE(1, h.max_chain_length());
            }

            EXPECT_FALSE( h.empty() );
            auto expected_size = max_ - (max_/2);
            EXPECT_EQ( expected_size, h.size());
            
            auto count = 0;
            h.forEach([&count](Value element){ ++ count; });
            EXPECT_EQ( expected_size, count);
            
            for(int i=max_-1; i>=max_/2; i--) {
                EXPECT_LE(1, h.max_chain_length());
                ASSERT_TRUE(h.has(i));
                h.remove(i);
                ASSERT_FALSE(h.has(i));
                EXPECT_TRUE(h.isConsistent());
            }
            
            EXPECT_TRUE(h.empty());
            EXPECT_EQ( 0, h.size() );
            EXPECT_EQ( 0, h.max_chain_length());
        }
        
    };

    
    template< typename Value >
    void run() {
        cout << "---" << endl;
        cout << "HashTable< "; COUT_TYPE(Value); cout << " >" << endl;
        cout << "---" << endl;
        
        Test<Value> test;

        test.run_logic();
    }
    
} // NS hashtable
} // NS test
} // NS imj


TEST(DataStructure, HashTable) {
    using namespace imj::test::hastable;

    run<int>();
}


