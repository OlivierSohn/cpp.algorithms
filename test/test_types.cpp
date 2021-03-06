
using namespace imajuscule;

TEST(Enum, initialization) {
    enum TestEnum {
        Test_v1=1,
        Test_v2,
    };
    
    TestEnum t{};
    EXPECT_EQ(0, static_cast<int>(t));
}

TEST(Power_of_2, funcs) {
    EXPECT_EQ(1, pow2(0));
    EXPECT_EQ(2, pow2(1));
    EXPECT_EQ(4, pow2(2));
    EXPECT_EQ(8, pow2(3));
    EXPECT_EQ(1024, pow2(10));

    EXPECT_EQ(0, power_of_two_exponent(0));
    EXPECT_EQ(0, power_of_two_exponent(1));
    EXPECT_EQ(1, power_of_two_exponent(2));
    EXPECT_EQ(2, power_of_two_exponent(4));
    EXPECT_EQ(3, power_of_two_exponent(8));
    
    EXPECT_EQ(0, ceil_power_of_two(0));
    EXPECT_EQ(1, ceil_power_of_two(1));
    EXPECT_EQ(2, ceil_power_of_two(2));
    EXPECT_EQ(4, ceil_power_of_two(3));
    EXPECT_EQ(4, ceil_power_of_two(4));
    EXPECT_EQ(8, ceil_power_of_two(5));
    EXPECT_EQ(8, ceil_power_of_two(6));
    EXPECT_EQ(8, ceil_power_of_two(7));
    EXPECT_EQ(8, ceil_power_of_two(8));
    EXPECT_EQ(16, ceil_power_of_two(9));

    EXPECT_EQ(0, floor_power_of_two(0));
    EXPECT_EQ(1, floor_power_of_two(1));
    EXPECT_EQ(2, floor_power_of_two(2));
    EXPECT_EQ(2, floor_power_of_two(3));
    EXPECT_EQ(4, floor_power_of_two(4));
    EXPECT_EQ(4, floor_power_of_two(5));
    EXPECT_EQ(4, floor_power_of_two(6));
    EXPECT_EQ(4, floor_power_of_two(7));
    EXPECT_EQ(8, floor_power_of_two(8));
    EXPECT_EQ(8, floor_power_of_two(9));
}

TEST(RelevantBits, funcs) {
    EXPECT_EQ(0, relevantBits(0));
    EXPECT_EQ(1, relevantBits(1));
    EXPECT_EQ(2, relevantBits(2));
    EXPECT_EQ(2, relevantBits(3));
    EXPECT_EQ(3, relevantBits(4));
    EXPECT_EQ(3, relevantBits(5));
}

TEST(Zeroes, LeadTrail) {
    EXPECT_EQ(0, count_trailing_zeroes(1));
    EXPECT_EQ(1, count_trailing_zeroes(2));
    EXPECT_EQ(31, count_leading_zeroes(1));
    EXPECT_EQ(30, count_leading_zeroes(2));
}

TEST(Cast, int) {
    ASSERT_EQ( std::numeric_limits<int>::max(), static_cast<int>(std::numeric_limits<int>::max()));
    auto i = static_cast<int>(.5f + std::numeric_limits<int>::max()/2);
    ASSERT_TRUE( i > 0 );
}
