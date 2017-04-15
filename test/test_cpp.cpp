
namespace imajuscule {
    namespace test_universal_ref {
        template<typename I>
        static void increment(I&& i) {
            ++i;
        }
    }
}
TEST(UniversalReference, pod) {
    int i=0;
    imajuscule::test_universal_ref::increment(i);
    ASSERT_EQ(1, i);
}


// reason for DISABLED :
// this test would pass if NaN values were "ignored" by std::minmax_element
TEST(MinMax, DISABLED_withNaN)
{
    
    {
        std::vector<float> v {0.f, 1.f, std::nanf("")};
        ASSERT_TRUE(v[2] != v[2]);
        auto res = std::minmax_element(v.begin(), v.end());
        ASSERT_EQ(0.f, *res.first);
        ASSERT_EQ(1.f, *res.second);
    }
    {
        std::vector<float> v {0.f, std::nanf(""), 1.f};
        ASSERT_TRUE(v[1] != v[1]);
        auto res = std::minmax_element(v.begin(), v.end());
        ASSERT_EQ(0.f, *res.first);
        ASSERT_EQ(1.f, *res.second);
    }
    {
        std::vector<float> v {std::nanf(""), 0.f, 1.f};
        ASSERT_TRUE(v[0] != v[0]);
        auto res = std::minmax_element(v.begin(), v.end());
        ASSERT_EQ(0.f, *res.first);
        ASSERT_EQ(1.f, *res.second);
    }
}
