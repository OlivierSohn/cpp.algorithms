
using namespace imajuscule;

TEST(Enum, initialization) {
    enum TestEnum {
        Test_v1=1,
        Test_v2,
    };
    
    TestEnum t{};
    EXPECT_EQ(0, static_cast<int>(t));
}
