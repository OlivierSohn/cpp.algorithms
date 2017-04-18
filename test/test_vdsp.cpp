
TEST(Accelerate, vmul) {
    float a[] {1.f, 2.f, 3.f, 4.f};
    float b[] {1.f, 2.f, 3.f, 4.f};
    float c[4];
    vDSP_vmul(a,1,b,1,c,1,4);
    ASSERT_FLOAT_EQ(1.f, c[0]);
    ASSERT_FLOAT_EQ(4.f, c[1]);
    ASSERT_FLOAT_EQ(9.f, c[2]);
    ASSERT_FLOAT_EQ(16.f, c[3]);
}

TEST(Accelerate, vmul_overlapp1) {
    float a[] {1.f, 2.f, 3.f, 4.f};
    float b[] {1.f, 2.f, 3.f, 4.f};
    vDSP_vmul(a,1,b,1,b,1,4);
    ASSERT_FLOAT_EQ(1.f, b[0]);
    ASSERT_FLOAT_EQ(4.f, b[1]);
    ASSERT_FLOAT_EQ(9.f, b[2]);
    ASSERT_FLOAT_EQ(16.f, b[3]);
}
TEST(Accelerate, vmul_overlapp2) {
    float a[] {1.f, 2.f, 3.f, 4.f};
    float b[] {1.f, 2.f, 3.f, 4.f};
    vDSP_vmul(a,1,b,1,a,1,4);
    ASSERT_FLOAT_EQ(1.f, a[0]);
    ASSERT_FLOAT_EQ(4.f, a[1]);
    ASSERT_FLOAT_EQ(9.f, a[2]);
    ASSERT_FLOAT_EQ(16.f, a[3]);
}

