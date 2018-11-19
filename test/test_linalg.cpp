
TEST(LinAlg, identity) {
  auto id = mkIdentity<float>(3);

  EXPECT_EQ(3, id[0].size());

  EXPECT_EQ(1, id[0][0]);
  EXPECT_EQ(1, id[1][1]);
  EXPECT_EQ(1, id[2][2]);

  EXPECT_EQ(0, id[0][1]);
  EXPECT_EQ(0, id[0][2]);
  EXPECT_EQ(0, id[1][0]);
  EXPECT_EQ(0, id[1][2]);
  EXPECT_EQ(0, id[2][0]);
  EXPECT_EQ(0, id[2][1]);
}

