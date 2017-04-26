

namespace testGradientDescent {
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
            auto index_min = 0;
            for(int i=1; i<values.size(); ++i) {
                if(values[i] < values[index_min]) {
                    index_min = i;
                }
            }
            
            int counter;
            auto f = create_f(values, counter);
            
            for(int i=0; i<values.size(); ++i) {
                counter = 0;
                float min_;
                int m = findLocalMinimum(n_iterations, i, f, min_);
                ASSERT_EQ(values[m], min_);
                ASSERT_EQ(values[index_min], values[m]);
#if DEBUG_GRADIENT_DESCENT == 0
                ASSERT_TRUE(counter <= 3 + std::abs(i - m)); // white box test
#endif
            }
            
            std::vector<int> invalid_starts {{
                -1, -2, -3, values.size(), values.size() + 1, values.size() + 2
            }};
            for(int i : invalid_starts) {
                float min_;
                ASSERT_THROW(findLocalMinimum(n_iterations, i, f, min_), std::logic_error) << i;
            }
        }
    }
}

TEST(GradientDescent, 1D) {
    using namespace testGradientDescent;
    
    // the inputs need to be concave for the test to work!
    
    // with unique global minimum:
    test({{
        0,1,2,3,4,5
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
    
    // with multiple global minimas:
    test({{
        1,1,2,3,4,5
    }});
    test({{
        1,1,1,3,4,5
    }});
    test({{
        1,1,1,1,4,5
    }});
    test({{
        1,1,1,1,1,5
    }});
    
    test({{
        1,1,1,1,1,1
    }});
    
    test({{
        5,4,3,2,1,0
    }});
    
    test({{
        5,4,3,2,0,0
    }});
    
    test({{
        5,4,3,0,0,0
    }});
    
    test({{
        5,4,0,0,0,0
    }});
    
    test({{
        5,0,0,0,0,0
    }});
    
    test({{
        0,0,0,0,0,0
    }});

}

