
TEST(Peaks, first_relevant_value) {
    using namespace imajuscule;
    std::vector<float> v{ 1.f, 2.f, 3.f };
    {
        auto it = first_relevant_value(v.begin(), v.end(), 3.5f);
        ASSERT_EQ(v.end(), it);
    }
    {
        auto it = first_relevant_value(v.begin(), v.end(), 2.5f);
        ASSERT_EQ(3.f, *it);
    }
    {
        auto it = first_relevant_value(v.begin(), v.end(), 1.5f);
        ASSERT_EQ(2.f, *it);
    }
    {
        auto it = first_relevant_value(v.begin(), v.end(), 0.5f);
        ASSERT_EQ(1.f, *it);
    }
    {
        auto it = first_relevant_value(v.begin(), v.end(), -0.5f);
        ASSERT_EQ(1.f, *it);
    }
}

TEST(Peaks, first_zero_crossing_forward)
{
    using namespace imajuscule;
    {
        std::vector<float> v{ 1.f, 2.f, 1.f };
        auto it = first_zero_crossing(v.begin(), v.end());
        ASSERT_EQ(v.end(), it);
    }
    {
        std::vector<float> v{ 1.f, 2.f, 1.f, 0.f };
        auto it = first_zero_crossing(v.begin(), v.end());
        ASSERT_EQ(0.f, *it);
    }
    {
        std::vector<float> v{ 1.f, 2.f, 0.f, 1.f };
        auto it = first_zero_crossing(v.begin(), v.end());
        ASSERT_EQ(0.f, *it);
    }
    {
        std::vector<float> v{ 1.f, 2.f, 0.f, -1.f };
        auto it = first_zero_crossing(v.begin(), v.end());
        ASSERT_EQ(0.f, *it);
    }
    {
        std::vector<float> v{ 1.f, 2.f, -1.f };
        auto it = first_zero_crossing(v.begin(), v.end());
        ASSERT_EQ(-1.f, *it);
    }
    {
        std::vector<float> v{ 1.f, 2.f, -1.f, 0.f };
        auto it = first_zero_crossing(v.begin(), v.end());
        ASSERT_EQ(-1.f, *it);
    }
}

TEST(Peaks, find_relevant_start)
{
    using namespace imajuscule;
    {
        std::vector<float> v{ -0.04f, -0.03f, -0.02f, -0.01f, 0.1f, 0.2f, 0.3f };
        auto it = find_relevant_start(v.begin(), v.end(), 0.15f);
        ASSERT_EQ(0.1f, *it);
    }
    {
        std::vector<float> v{ -0.04f, -0.03f, -0.02f, -0.01f, 0.1f, 0.2f, 0.3f };
        auto it = find_relevant_start(v.begin(), v.end(), 0.25f);
        // the algorithm finds the first relevant value, and goes backward from there to find the first sign change and returns the sample just after
        ASSERT_EQ(0.1f, *it);
    }
    {
        std::vector<float> v{ 0.04f, 0.03f, 0.02f, 0.01f, -0.1f, -0.2f, -0.3f };
        auto it = find_relevant_start(v.begin(), v.end(), 0.25f);
        // the algorithm finds the first relevant value, and goes backward from there to find the first sign change and returns the sample just after
        ASSERT_EQ(-0.1f, *it);
    }
    {
        std::vector<float> v{ -0.04f, -0.03f, -0.02f, -0.01f, 0.1f, 0.2f, 0.3f };
        auto it = find_relevant_start(v.begin(), v.end(), 0.5f);
        // the level is too high and was never reached so "end" should be returned
        ASSERT_EQ(v.end(), it);
    }
    {
        std::vector<float> v{ 1.f, 2.f, 1.f };
        auto it = find_relevant_start(v.begin(), v.end(), 0.5f);
        ASSERT_EQ(v.begin(), it);
    }
    {
        std::vector<float> v{ -1.f, -2.f, -1.f };
        auto it = find_relevant_start(v.begin(), v.end(), 0.5f);
        ASSERT_EQ(v.begin(), it);
    }
    {
        std::vector<float> v;
        auto it = find_relevant_start(v.begin(), v.end(), 0.5f);
        ASSERT_EQ(v.end(), it);
    }
    {
        std::vector<float> v{.1f};
        auto it = find_relevant_start(v.begin(), v.end(), 0.5f);
        ASSERT_EQ(v.end(), it);
    }
    {
        std::vector<float> v{1.f};
        auto it = find_relevant_start(v.begin(), v.end(), 0.5f);
        ASSERT_EQ(1.f, *it);
    }
    {
        std::vector<float> v{-1.f};
        auto it = find_relevant_start(v.begin(), v.end(), 0.5f);
        ASSERT_EQ(-1.f, *it);
    }
}

TEST(Peaks, first_non_abs_decreasing)
{
    using namespace imajuscule;
    {
        std::vector<float> v{ 0.f, 1.f };
        auto it = first_non_abs_decreasing(v.begin(), v.end());
        ASSERT_EQ(0.f, *it);
    }
    {
        std::vector<float> v{ 1.f, 0.f };
        auto it = first_non_abs_decreasing(v.begin(), v.end());
        ASSERT_EQ(0.f, *it);
    }
    {
        std::vector<float> v{ 0.f, -1.f };
        auto it = first_non_abs_decreasing(v.begin(), v.end());
        ASSERT_EQ(0.f, *it);
    }
    {
        std::vector<float> v{ -1.f, 0.f };
        auto it = first_non_abs_decreasing(v.begin(), v.end());
        ASSERT_EQ(0.f, *it);
    }
    {
        std::vector<float> v{ 3.f, 2.f, 1.f, 2.f, 3.f };
        auto it = first_non_abs_decreasing(v.begin(), v.end());
        ASSERT_EQ(1.f, *it);
    }
    {
        std::vector<float> v{ -3.f, -2.f, -1.f, -2.f, -3.f };
        auto it = first_non_abs_decreasing(v.begin(), v.end());
        ASSERT_EQ(-1.f, *it);
    }
}


TEST(Peaks, first_non_abs_avg_decreasing_unit)
{
    using namespace imajuscule;
    constexpr auto sliding_avg_size = 1;
    {
        std::vector<float> v{ 0.f, 1.f };
        auto it = first_non_abs_avg_decreasing(v.begin(), v.end(), sliding_avg_size);
        ASSERT_EQ(0.f, *it);
    }
    {
        std::vector<float> v{ 1.f, 0.f };
        auto it = first_non_abs_avg_decreasing(v.begin(), v.end(), sliding_avg_size);
        ASSERT_EQ(0.f, *it);
    }
    {
        std::vector<float> v{ 0.f, -1.f };
        auto it = first_non_abs_avg_decreasing(v.begin(), v.end(), sliding_avg_size);
        ASSERT_EQ(0.f, *it);
    }
    {
        std::vector<float> v{ -1.f, 0.f };
        auto it = first_non_abs_avg_decreasing(v.begin(), v.end(), sliding_avg_size);
        ASSERT_EQ(0.f, *it);
    }
    {
        std::vector<float> v{ 3.f, 2.f, 1.f, 2.f, 3.f };
        auto it = first_non_abs_avg_decreasing(v.begin(), v.end(), sliding_avg_size);
        ASSERT_EQ(1.f, *it);
    }
    {
        std::vector<float> v{ -3.f, -2.f, -1.f, -2.f, -3.f };
        auto it = first_non_abs_avg_decreasing(v.begin(), v.end(), sliding_avg_size);
        ASSERT_EQ(-1.f, *it);
    }
}

TEST(Peaks, first_non_abs_avg_decreasing)
{
    using namespace imajuscule;
    constexpr auto sliding_avg_size = 2;
    {
        std::vector<float> v{ 0.f, 1.f };
        auto it = first_non_abs_avg_decreasing(v.begin(), v.end(), sliding_avg_size);
        ASSERT_EQ(0.f, *it);
    }
    {
        std::vector<float> v{ 1.f, 0.f };
        auto it = first_non_abs_avg_decreasing(v.begin(), v.end(), sliding_avg_size);
        ASSERT_EQ(0.f, *it);
    }
    {
        std::vector<float> v{ 0.f, -1.f };
        auto it = first_non_abs_avg_decreasing(v.begin(), v.end(), sliding_avg_size);
        ASSERT_EQ(0.f, *it);
    }
    {
        std::vector<float> v{ -1.f, 0.f };
        auto it = first_non_abs_avg_decreasing(v.begin(), v.end(), sliding_avg_size);
        ASSERT_EQ(0.f, *it);
    }
    {
        std::vector<float> v{ 3.f, 2.f, 1.f, 1.5f, 2.f, 2.5f, 3.f };
        auto it = first_non_abs_avg_decreasing(v.begin(), v.end(), sliding_avg_size);
        ASSERT_EQ(1.5f, *it);
    }
    {
        std::vector<float> v{ -3.f, -2.f, -1.f, -1.5f, -2.f, -2.5f, -3.f };
        auto it = first_non_abs_avg_decreasing(v.begin(), v.end(), sliding_avg_size);
        ASSERT_EQ(-1.5f, *it);
    }
}

TEST(Peaks, max_abs_integrated_lobe)
{
    using namespace imajuscule;
    {
        std::vector<float> v{ 0.f, 1.f, 0.f };
        auto val = max_abs_integrated_lobe(v.begin(), v.end());
        ASSERT_EQ(1.f, val);
    }
    {
        std::vector<float> v{ 0.f, -1.f, 0.f };
        auto val = max_abs_integrated_lobe(v.begin(), v.end());
        ASSERT_EQ(1.f, val);
    }
    {
        std::vector<float> v{ 1.f };
        auto val = max_abs_integrated_lobe(v.begin(), v.end());
        ASSERT_EQ(1.f, val);
    }
    {
        std::vector<float> v{ -1.f };
        auto val = max_abs_integrated_lobe(v.begin(), v.end());
        ASSERT_EQ(1.f, val);
    }
    {
        std::vector<float> v{ 0.f, 1.f, 2.f, 0.f };
        auto val = max_abs_integrated_lobe(v.begin(), v.end());
        ASSERT_EQ(3.f, val);
    }
    {
        std::vector<float> v{ 0.f, 1.f, 2.f, -4.f };
        auto val = max_abs_integrated_lobe(v.begin(), v.end());
        ASSERT_EQ(4.f, val);
    }
    {
        std::vector<float> v{ -4.f, 0.f, 1.f, 2.f };
        auto val = max_abs_integrated_lobe(v.begin(), v.end());
        ASSERT_EQ(4.f, val);
    }
    {
        std::vector<float> v{ 3.f, -4.f, 1.f, 2.f, -3.f };
        auto val = max_abs_integrated_lobe(v.begin(), v.end());
        ASSERT_EQ(4.f, val);
    }
}

TEST(Peaks, abs_integrated)
{
    using namespace imajuscule;
    {
        std::vector<float> v{ 0.f, 1.f, 0.f };
        auto val = abs_integrated(v.begin(), v.end());
        ASSERT_EQ(1.f, val);
    }
    {
        std::vector<float> v{ 0.f, -1.f, 0.f };
        auto val = abs_integrated(v.begin(), v.end());
        ASSERT_EQ(1.f, val);
    }
    {
        std::vector<float> v{ 0.f, -1.f, 1.f };
        auto val = abs_integrated(v.begin(), v.end());
        ASSERT_EQ(2.f, val);
    }
}
