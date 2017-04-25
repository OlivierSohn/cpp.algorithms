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

TEST(Math, expMean) {
    ASSERT_EQ(0, exp_mean(0,0));
    ASSERT_EQ(3, exp_mean(3,3));
    ASSERT_EQ(-3, exp_mean(-3,-3));
    
    auto v = exp_mean(2,128);
    ASSERT_EQ(16, v);
}
