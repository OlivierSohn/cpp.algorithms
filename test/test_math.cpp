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

TEST(PhasedSum, test) {
    using namespace imajuscule;
    constexpr auto cyclic_size = 5;
    cyclic<float> c(cyclic_size);
    c.feed(1);
    c.feed(2);
    c.feed(3);
    c.feed(4);
    c.feed(5);

    cyclic<float> res(cyclic_size);
    {
        constexpr auto phase = 0;
        constexpr auto n_iterators = 2;
        compute_phased_sum(c, phase, n_iterators, res);

        auto it = res.begin();
        ASSERT_FLOAT_EQ(2, *it);
        ++it;
        ASSERT_FLOAT_EQ(4, *it);
        ++it;
        ASSERT_FLOAT_EQ(6, *it);
        ++it;
        ASSERT_FLOAT_EQ(8, *it);
        ++it;
        ASSERT_FLOAT_EQ(10, *it);
        ++it;
    }
    {
        constexpr auto phase = 1;
        constexpr auto n_iterators = 2;
        compute_phased_sum(c, phase, n_iterators, res);

        auto it = res.begin();
        ASSERT_FLOAT_EQ(3, *it);
        ++it;
        ASSERT_FLOAT_EQ(5, *it);
        ++it;
        ASSERT_FLOAT_EQ(7, *it);
        ++it;
        ASSERT_FLOAT_EQ(9, *it);
        ++it;
        ASSERT_FLOAT_EQ(6, *it);
        ++it;
    }

    for(int i=0; i<10; ++i) {
        auto phase = 1 + cyclic_size;
        constexpr auto n_iterators = 2;
        compute_phased_sum(c, phase, n_iterators, res);

        auto it = res.begin();
        ASSERT_FLOAT_EQ(3, *it);
        ++it;
        ASSERT_FLOAT_EQ(5, *it);
        ++it;
        ASSERT_FLOAT_EQ(7, *it);
        ++it;
        ASSERT_FLOAT_EQ(9, *it);
        ++it;
        ASSERT_FLOAT_EQ(6, *it);
        ++it;
    }

    {
        constexpr auto phase = 2;
        constexpr auto n_iterators = 3;
        compute_phased_sum(c, phase, n_iterators, res);

        auto it = res.begin();
        ASSERT_FLOAT_EQ(9, *it);
        ++it;
        ASSERT_FLOAT_EQ(7, *it);
        ++it;
        ASSERT_FLOAT_EQ(10, *it);
        ++it;
        ASSERT_FLOAT_EQ(8, *it);
        ++it;
        ASSERT_FLOAT_EQ(11, *it);
        ++it;
    }

    for(int i=0; i<10; ++i) {
        auto phase = 2 + cyclic_size;
        constexpr auto n_iterators = 3;
        compute_phased_sum(c, phase, n_iterators, res);

        auto it = res.begin();
        ASSERT_FLOAT_EQ(9, *it);
        ++it;
        ASSERT_FLOAT_EQ(7, *it);
        ++it;
        ASSERT_FLOAT_EQ(10, *it);
        ++it;
        ASSERT_FLOAT_EQ(8, *it);
        ++it;
        ASSERT_FLOAT_EQ(11, *it);
        ++it;
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

TEST(Math, IntDiv) {
    for(int i=1; i<10; ++i) {
        for(int j=1; j<10; ++j) {
            ASSERT_EQ(j, (j/i)*i + (j%i));
        }
    }
}


TEST(UnsignedCyclicDist, simple)
{
  using namespace imajuscule;
  {
    unsigned int a = 3;
    unsigned int b = 5;
    ASSERT_EQ(cyclic_unsigned_dist(a,b), 2);
    ASSERT_EQ(cyclic_unsigned_dist(b,a), 2);
  }
  {
    unsigned int a = std::numeric_limits<unsigned int>::min();
    unsigned int b = std::numeric_limits<unsigned int>::max();
    ASSERT_EQ(cyclic_unsigned_dist(a,b), 1);
    ASSERT_EQ(cyclic_unsigned_dist(b,a), 1);
  }
  {
    unsigned int a = std::numeric_limits<unsigned int>::min() + 2;
    unsigned int b = std::numeric_limits<unsigned int>::max() - 3;
    ASSERT_EQ(cyclic_unsigned_dist(a,b), 6);
    ASSERT_EQ(cyclic_unsigned_dist(b,a), 6);
  }
  unsigned int const one = 1;
  unsigned int const zero = 0;
  ASSERT_EQ(cyclic_unsigned_dist(one,one), 0);
  ASSERT_EQ(cyclic_unsigned_dist(zero,zero), 0);
  ASSERT_EQ(cyclic_unsigned_dist(std::numeric_limits<unsigned int>::max(),std::numeric_limits<unsigned int>::max()), 0);
}
