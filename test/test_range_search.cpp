
namespace testRangeSearch {
    template<typename ARRAY>
    auto create_f(ARRAY const & array, int & counter) {
        return [&array, &counter](int x, float & val) {
            using namespace imajuscule;
            ++counter;
            if(x < 0 || x >= array.size()) {
                return ParamState::OutOfRange;
            }
            val = array[x];
            return ParamState::Ok;
        };
    }
    
    void test(std::array<int, 6> values) {
        using namespace imajuscule;
        
        auto n_iterations = 1;
        for(; n_iterations < 4; ++n_iterations) {
            int counter;
            auto f = create_f(values, counter);
            
            for(int i=0; i<values.size(); ++i) {
                counter = 0;
                float min_;
                int m = findGlobalMinimum({ 0, 5 },
                                          f,
                                          min_);
                ASSERT_EQ(0, min_);
                ASSERT_EQ(values[m], min_);
                ASSERT_EQ(6, counter);
            }
            
            for(int i=0; i<values.size(); ++i) {
                counter = 0;
                float min_;
                constexpr auto threshold = 8;
                int m = findGlobalMinimumOrSmallerThan({ 0, 5 },
                                                       threshold,
                                                       f,
                                                       min_);
                ASSERT_EQ(values[m], min_);
                ASSERT_TRUE(values[m] < threshold);
                ASSERT_EQ(1, counter);
            }

            
            for(int i=0; i<values.size(); ++i) {
                counter = 0;
                float min_;
                constexpr auto threshold = 3;
                int m = findGlobalMinimumOrSmallerThan({ 0, 5 },
                                                       threshold,
                                                       f,
                                                       min_);
                ASSERT_EQ(values[m], min_);
                ASSERT_TRUE(values[m] < threshold);
                ASSERT_TRUE(counter <= 4);
            }
        }
    }
}

TEST(RangeSearch, 1D) {
    using namespace testRangeSearch;
    
    test({{
        0,1,2,3,4,5
    }});
    
    test({{
        1,0,5,3,4,2
    }});
    
    test({{
        5,4,3,2,1,0
    }});
    
    test({{
        5,0,1,2,3,4
    }});
    
    test({{
        5,4,0,1,2,3
    }});
    
    test({{
        5,4,3,0,1,2
    }});
    test({{
        5,4,3,2,0,1
    }});
    test({{
        5,4,3,2,1,0
    }});
}

