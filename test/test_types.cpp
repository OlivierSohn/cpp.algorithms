
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
    EXPECT_EQ(0, power_of_two_exponent(1));
    EXPECT_EQ(1, power_of_two_exponent(2));
    EXPECT_EQ(2, power_of_two_exponent(4));
    EXPECT_EQ(3, power_of_two_exponent(8));
}

TEST(RelevantBits, funcs) {
    EXPECT_EQ(0, relevantBits(0));
    EXPECT_EQ(1, relevantBits(1));
    EXPECT_EQ(2, relevantBits(2));
    EXPECT_EQ(2, relevantBits(3));
    EXPECT_EQ(3, relevantBits(4));
    EXPECT_EQ(3, relevantBits(5));
}

