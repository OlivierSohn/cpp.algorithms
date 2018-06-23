
TEST(PairArrayDistant, simple) {
  using namespace imajuscule;

  DistantPairArray<int,float> pa(10,2,4.f);

  EXPECT_EQ(10, pa.size());

  for(int i=0; i<10; ++i) {
    EXPECT_EQ(2, pa.firsts()[i]);
    EXPECT_FLOAT_EQ(4.f, pa.seconds()[i]);
  }
}


TEST(PairArrayDistant, AsBs) {
  using namespace imajuscule;
  DistantPairArray<BigAlignment,SmallAlignment> pa(10, BigAlignment{1.0}, SmallAlignment{2.0});
  static_assert(decltype(pa)::order == detail::Order::As_Bs);


  EXPECT_EQ(10, pa.size());
  
  for(int i=0; i<10; ++i) {
    EXPECT_FLOAT_EQ(1.0, pa.firsts()[i].v);
    EXPECT_FLOAT_EQ(2.0, pa.seconds()[i].v);
  }
}

TEST(PairArrayDistant, BsAs) {
  using namespace imajuscule;
  DistantPairArray<SmallAlignment,BigAlignment> pa(10, SmallAlignment{1.0}, BigAlignment{2.0});
  static_assert(decltype(pa)::order == detail::Order::Bs_As);

  EXPECT_EQ(10, pa.size());
  
  for(int i=0; i<10; ++i) {
    EXPECT_FLOAT_EQ(1.0, pa.firsts()[i].v);
    EXPECT_FLOAT_EQ(2.0, pa.seconds()[i].v);
  }
}

TEST(PairArrayDistant, RangeIteration) {
  using namespace imajuscule;
  DistantPairArray<int8_t,std::string> pa(10, 1, "hello");

  {
    int n = 0;
    for(auto & p : firsts(pa)) {
      EXPECT_EQ(1, p);
      ++n;
    }
    EXPECT_EQ(10, n);
  }
  {
    int n = 0;
    for(auto & p : seconds(pa)) {
      EXPECT_EQ("hello", p);
      ++n;
    }
    EXPECT_EQ(10, n);
  }
}

TEST(PairArrayDistant, RangeIteration2) {
  using namespace imajuscule;
  DistantPairArray<std::string,int8_t> pa(10, "hello", 1);
  
  for(auto & p : firsts(pa)) {
    EXPECT_EQ("hello", p);
  }
  for(auto & p : seconds(pa)) {
    EXPECT_EQ(1, p);
  }
}
