TEST(Math, dichotomic_sum) {
    using namespace imajuscule;
    
    {
        std::vector<float> v{};
        ASSERT_NEAR(0, dichotomic_sum(v.begin(), v.end()), 1e-5);
    }
    
    {
        std::vector<float> v{0};
        ASSERT_NEAR(0, dichotomic_sum(v.begin(), v.end()), 1e-5);
    }
    
    {
        std::vector<float> v{1};
        ASSERT_NEAR(1, dichotomic_sum(v.begin(), v.end()), 1e-5);
    }
    
    {
        std::vector<float> v{{1,2}};
        ASSERT_NEAR(3, dichotomic_sum(v.begin(), v.end()), 1e-5);
    }
    
    {
        std::vector<float> v{{2,3,4}};
        ASSERT_NEAR(9, dichotomic_sum(v.begin(), v.end()), 1e-5);
    }
}

TEST(MaxSlidingSum, test) {
    using namespace imajuscule;
    constexpr auto cyclic_size = 5;
    cyclic<float> c(cyclic_size);
    c.feed(1);
    c.feed(5);
    c.feed(2);
    c.feed(4);
    c.feed(3);
    
    float full_sum = 1+2+3+4+5;
    float max_value = 5;
    
    std::array<float,cyclic_size + 1> expected_max{{
        0,
        5,
        5 + 2,
        5 + 2 + 4,
        5 + 2 + 4 + 3,
        5 + 2 + 4 + 3 + 1
    }};
    
    // test with slidingWidth == a number of periods of the cycle
    for(int p=1; p<4; ++p) {
        ASSERT_FLOAT_EQ(p * full_sum,
                        computeMaxSlidingSum(c, p*cyclic_size));
        
    }
    
    // test with slidingWidth <=  size of the cycle
    for(int i=1; i <= cyclic_size; ++i) {
        ASSERT_FLOAT_EQ(expected_max[i],
                        computeMaxSlidingSum(c, i)) << i;
    }
    
    // test with slidingWidth >  size of the cycle
    for(auto p = 1; p<4; ++p) {
        for(int i=0; i<=cyclic_size; ++i) {
            ASSERT_FLOAT_EQ(expected_max[i] + p*full_sum,
                            computeMaxSlidingSum(c, p*cyclic_size + i)) << p;
        }
    }
}

TEST(Math, expMean) {
    using namespace imajuscule;
    
    ASSERT_EQ(0, exp_mean(0,0));
    ASSERT_EQ(3, exp_mean(3,3));
    ASSERT_EQ(-3, exp_mean(-3,-3));
    
    auto v = exp_mean(2,128);
    ASSERT_EQ(16, v);
}
