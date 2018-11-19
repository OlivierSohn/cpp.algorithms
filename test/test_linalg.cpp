
TEST(LinAlg, identity) {
  using namespace imajuscule;
  
  auto id = mkIdentity<float>(3);

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

TEST(LinAlg, croutLU) {
  using namespace imajuscule;
  Matrix<float> A, L, U;
  constexpr int sz = 4;
  A.resize(sz, sz);
  L.resize(sz, sz);
  U.resize(sz, sz);

  // using random numbers so that A is invertible
  for(int i=0; i<sz; ++i) {
    for(int j=0; j<sz; ++j) {
      A[i][j] = std::uniform_real_distribution<>{
        0.f, 1.f
      }(mersenne<SEEDED::No>());
    }
  }
  
  // crout performs the LU decomposition of A
  auto res = crout(A,L,U);
  
  ASSERT_TRUE(res); // else A may not be invertible
  
  ASSERT_TRUE(isUpper(U));
  ASSERT_TRUE(isLower(L));
  
  // verify that 'A = LU'
  {
    Matrix<float> Compare;
    Compare.resize(sz, sz);
    multiply(L,U,Compare);
    ASSERT_TRUE(equals(A, Compare, 0.001f));
  }
}

TEST(LinAlg, solve) {
  using namespace imajuscule;

  // let's solve this:
  
  //  x + 2y - 3z = 1
  // 3x -  y +  z = 5
  // 5x + 3y - 2z = 7
  
  // x = -2y +3z + 1
  // -6y + 9z +3 -y + z = 5   ->   -7y +10z = 2   (a)
  // -10y +15z +5 +3y -2z = 7 ->   -7y +13z = 2  (b)

  // (a)-(b): -3z = 0 ->                   z = 0
  // -7y = 2 ->           y = -2/7
  // x = 11/7
  //
  // the system is equivalent to 'Av=w' where:
  //
  // w = T[1 5 7]
  //
  // v = T[x y z]
  //
  // A = [1 2 -3
  //      3 -1 1
  //      5 3 -2]

  Matrix<double> A;
  A.resize(3,3);
  A[0][0] = 1.; A[0][1] =  2.; A[0][2] = -3.;
  A[1][0] = 3.; A[1][1] = -1.; A[1][2] =  1.;
  A[2][0] = 5.; A[2][1] =  3.; A[2][2] = -2.;
  
  std::vector<double> w{1., 5., 7.};
  std::vector<double> v;

  auto res = solve(A, w, v);
  
  ASSERT_TRUE(res);
  ASSERT_EQ(3, v.size());
  ASSERT_EQ(11./7., v[0]);
  ASSERT_EQ(-2./7., v[1]);
  ASSERT_EQ(0., v[2]);
}
