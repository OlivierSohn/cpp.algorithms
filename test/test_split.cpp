TEST(split, simple) {
  using namespace imajuscule;
  using S = Split<uint64_t, 54, 10>;
  static_assert(S::maxLow  == (pow2(10) - 1));
  static_assert(S::maxHigh == (pow2(54) - 1));
  {
    S split(std::numeric_limits<uint64_t>::max());
    EXPECT_EQ(S::maxLow, split.getLow());
    EXPECT_EQ(S::maxHigh, split.getHigh());
    EXPECT_EQ(S::maxHigh << 10, split.getHighWithZeros());
  }
  {
    S split(0);
    EXPECT_EQ(0, split.getLow());
    EXPECT_EQ(0, split.getHigh());
    EXPECT_EQ(0, split.getHighWithZeros());
  }
  std::vector<uint64_t> someLows, someHighs;
  for(uint64_t i=0; i<=5; ++i) {
    someLows.push_back(i);
    someLows.push_back(S::maxLow-i);
  }
  for(uint64_t i=0; i<=5; ++i) {
    someHighs.push_back(i);
    someHighs.push_back(S::maxHigh - i);
  }
  
  for(uint64_t i : someLows)
  {
    S split(i);
    EXPECT_EQ(i, split.getLow());
    EXPECT_EQ(0, split.getHigh());
    EXPECT_EQ(0, split.getHighWithZeros());
  }
  for(uint64_t i : someHighs)
  {
    S split(i << 10);
    EXPECT_EQ(0, split.getLow());
    EXPECT_EQ(i, split.getHigh());
    EXPECT_EQ(i << 10, split.getHighWithZeros());
  }

  for(uint64_t low : someLows)
  {
    for(uint64_t high : someHighs)
    {
      S split((high << 10) + low);
      EXPECT_EQ(low, split.getLow());
      EXPECT_EQ(high, split.getHigh());
      EXPECT_EQ(high << 10, split.getHighWithZeros());
    }
  }
}

TEST(split, complex) {
  using namespace imajuscule;
  using S = Split<uint64_t, 54, 10>;
  {
    // verify asserts are thrown when the input arguments are out of bound:
    auto f = [](auto a, auto b){
      S tmp(a,b);
    };
    auto g = [=]() {
      f(0,S::maxLow+1);
    };
    auto h = [=]() {
      f(S::maxHigh+1,0);
    };

    EXPECT_DEBUG_ASSERT(g());
    EXPECT_DEBUG_ASSERT(h());
  }
  {
    S split(S::maxHigh, S::maxLow);
    EXPECT_EQ(S::maxLow, split.getLow());
    EXPECT_EQ(S::maxHigh, split.getHigh());
    EXPECT_EQ(S::maxHigh << 10, split.getHighWithZeros());
  }
  {
    S split(0,0);
    EXPECT_EQ(0, split.getLow());
    EXPECT_EQ(0, split.getHigh());
    EXPECT_EQ(0, split.getHighWithZeros());
  }
  std::vector<uint64_t> someLows, someHighs;
  for(uint64_t i=0; i<=5; ++i) {
    someLows.push_back(i);
    someLows.push_back(S::maxLow-i);
  }
  for(uint64_t i=0; i<=5; ++i) {
    someHighs.push_back(i);
    someHighs.push_back(S::maxHigh - i);
  }
  
  for(uint64_t i : someLows)
  {
    S split(0,i);
    EXPECT_EQ(i, split.getLow());
    EXPECT_EQ(0, split.getHigh());
    EXPECT_EQ(0, split.getHighWithZeros());
  }
  for(uint64_t i : someHighs)
  {
    S split(i,0);
    EXPECT_EQ(0, split.getLow());
    EXPECT_EQ(i, split.getHigh());
    EXPECT_EQ(i << 10, split.getHighWithZeros());
  }
  
  for(uint64_t low : someLows)
  {
    for(uint64_t high : someHighs)
    {
      S split(high,low);
      EXPECT_EQ(low, split.getLow());
      EXPECT_EQ(high, split.getHigh());
      EXPECT_EQ(high << 10, split.getHighWithZeros());
    }
  }
  
}
