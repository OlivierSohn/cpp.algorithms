

TEST(SlidingAverage, test) {
    using namespace imajuscule;
    {
        slidingAverage<float> avg(3);
        EXPECT_FLOAT_EQ(0.f, avg.compute());
        EXPECT_FLOAT_EQ(0.f, avg.compute());
        avg.feed(3.f);
        EXPECT_EQ(3, avg.size());
        EXPECT_FLOAT_EQ(3.f, avg.compute());
        EXPECT_FLOAT_EQ(3.f, avg.compute());
        avg.feed(3.f);
        EXPECT_EQ(3, avg.size());
        EXPECT_FLOAT_EQ(3.f, avg.compute());
        EXPECT_FLOAT_EQ(3.f, avg.compute());
        avg.feed(0.f);
        EXPECT_EQ(3, avg.size());
        EXPECT_FLOAT_EQ(2.f, avg.compute());
        EXPECT_FLOAT_EQ(2.f, avg.compute());
        avg.feed(0.f);
        EXPECT_EQ(3, avg.size());
        EXPECT_FLOAT_EQ(1.f, avg.compute());
        EXPECT_FLOAT_EQ(1.f, avg.compute());
        avg.feed(0.f);
        EXPECT_EQ(3, avg.size());
        EXPECT_FLOAT_EQ(0.f, avg.compute());
        EXPECT_FLOAT_EQ(0.f, avg.compute());
        avg.feed(0.f);
        EXPECT_EQ(3, avg.size());
        EXPECT_FLOAT_EQ(0.f, avg.compute());
        EXPECT_FLOAT_EQ(0.f, avg.compute());
    }
    {
        slidingAverage<float, CyclicInitialization::INITIAL_VALUES> avg(3);
        EXPECT_FLOAT_EQ(0.f, avg.compute());
        EXPECT_FLOAT_EQ(0.f, avg.compute());
        avg.feed(3.f);
        EXPECT_EQ(3, avg.size());
        EXPECT_FLOAT_EQ(1.f, avg.compute());
        EXPECT_FLOAT_EQ(1.f, avg.compute());
        avg.feed(3.f);
        EXPECT_EQ(3, avg.size());
        EXPECT_FLOAT_EQ(2.f, avg.compute());
        EXPECT_FLOAT_EQ(2.f, avg.compute());
        avg.feed(3.f);
        EXPECT_EQ(3, avg.size());
        EXPECT_FLOAT_EQ(3.f, avg.compute());
        EXPECT_FLOAT_EQ(3.f, avg.compute());
        avg.feed(3.f);
        EXPECT_EQ(3, avg.size());
        EXPECT_FLOAT_EQ(3.f, avg.compute());
        EXPECT_FLOAT_EQ(3.f, avg.compute());
        avg.feed(0.f);
        EXPECT_EQ(3, avg.size());
        EXPECT_FLOAT_EQ(2.f, avg.compute());
        EXPECT_FLOAT_EQ(2.f, avg.compute());
        avg.feed(0.f);
        EXPECT_EQ(3, avg.size());
        EXPECT_FLOAT_EQ(1.f, avg.compute());
        EXPECT_FLOAT_EQ(1.f, avg.compute());
        avg.feed(0.f);
        EXPECT_EQ(3, avg.size());
        EXPECT_FLOAT_EQ(0.f, avg.compute());
        EXPECT_FLOAT_EQ(0.f, avg.compute());
        avg.feed(0.f);
        EXPECT_EQ(3, avg.size());
        EXPECT_FLOAT_EQ(0.f, avg.compute());
        EXPECT_FLOAT_EQ(0.f, avg.compute());
    }
}


TEST(SlidingWindowAverage, test) {
    using namespace imajuscule;
    {
        slidingWindowedAverage<float> a(3, [](auto ratio){ return ratio; });
        a.feed(1);
        a.feed(1);
        a.feed(1);
        ASSERT_NEAR(2.f/3.f, a.compute(), 1e-6);
    }
    {
        slidingWindowedAverage<float> a(3, [](auto ratio){ return ratio; });
        a.feed(0);
        a.feed(1);
        a.feed(0);
        ASSERT_NEAR(1.f/3.f, a.compute(), 1e-6);
    }
    {
        slidingWindowedAverage<float> a(3, [](auto ratio){ return ratio; });
        a.feed(1);
        a.feed(0);
        a.feed(0);
        ASSERT_NEAR(.5f/3.f, a.compute(), 1e-6);
    }
    {
        slidingWindowedAverage<float> a(4, [](auto ratio){ return ratio; });
        a.feed(1);
        a.feed(1);
        a.feed(1);
        a.feed(1);
        ASSERT_NEAR(.6f, a.compute(), 1e-6);
    }
    {
        slidingWindowedAverage<float> a(4, [](auto ratio){ return ratio; });
        a.feed(0);
        a.feed(1);
        a.feed(1);
        a.feed(0);
        ASSERT_NEAR(.4f, a.compute(), 1e-6);
    }
    {
        slidingWindowedAverage<float> a(4, [](auto ratio){ return ratio; });
        a.feed(1);
        a.feed(0);
        a.feed(0);
        a.feed(1);
        ASSERT_NEAR(.2f, a.compute(), 1e-6);
    }
    {
        slidingWindowedAverage<float> a(4, [](auto ratio){ return ratio; });
        a.feed(0);
        a.feed(0);
        a.feed(0);
        a.feed(1);
        ASSERT_NEAR(.2f/2, a.compute(), 1e-6);
    }
    {
        slidingWindowedAverage<float> a(4, [](auto ratio){ return ratio; });
        a.feed(1);
        a.feed(0);
        a.feed(0);
        a.feed(0);
        ASSERT_NEAR(.2f/2, a.compute(), 1e-6);
    }
}
