
TEST(math,angles) {
    using namespace imajuscule;
    EXPECT_NEAR(0.f, modulo_angle(2*M_PI), 1e-6);
    EXPECT_NEAR(0.f, modulo_angle(4*M_PI), 1e-6);
    EXPECT_NEAR(0.f, modulo_angle(0.f), 1e-6);
    auto pi = modulo_angle(M_PI);
    if(pi > 0) {
        EXPECT_NEAR(M_PI, pi, 1e-6);
    }
    else {
        EXPECT_NEAR(-M_PI, pi, 1e-6);
    }
}
