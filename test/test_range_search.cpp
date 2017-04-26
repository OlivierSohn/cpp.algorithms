
namespace testRangeSearch {

    /*
     Creates a "V" function where the bottom of the V touches the x axis at 'x==a'
     */
    constexpr auto n_lin_steps = 256;
    auto create_lin_f(int min_step, int & counter) {
        auto a = min_step / static_cast<float>(n_lin_steps-1);
        return [a, &counter](int p_x, float & val) {
            assert(p_x >= 0);
            assert(p_x < n_lin_steps);
            
            using namespace imajuscule;
            ++counter;
            
            auto x = p_x / static_cast<float>(n_lin_steps-1);
            
            if(a < 0.00001) {
                val = x;
            }
            else if(a > 0.99999) {
                val = 1-x;
            }
            else if(x < a) {
                val = 1 - x/a;
            }
            else if(a != 1) {
                val = (x-a)/(1-a);
            }
            else {
                val = 1;
            }
            
            return ParamState::Ok;
        };
    }
    
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
        
        int counter;
        auto f = create_f(values, counter);
        
        for(int i=0; i<values.size(); ++i) {
            counter = 0;
            float min_;
            int m = findRangedLocalMinimum(
                                           1,
                                           { 0, 5 },
                                           f,
                                           min_);
            ASSERT_EQ(0, min_);
            ASSERT_EQ(values[m], min_);
            ASSERT_TRUE(6 == counter || 5 == counter);
        }
    }
    
    void test_lin() {
        using namespace imajuscule;
        
        int counter;
        for(auto m = 0; m < n_lin_steps; ++m) {
            auto f = create_lin_f(m, counter);
            
            counter = 0;
            float min_;
            int m_index = findRangedLocalMinimum(
                                                 1,
                                                 { 0, n_lin_steps-1 },
                                                 f,
                                                 min_);
            ASSERT_EQ(0.f, min_);
            ASSERT_EQ(m, m_index);
        }
    }
}

TEST(RangeSearch, 1D) {
    using namespace testRangeSearch;
    
    test_lin();
    
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

