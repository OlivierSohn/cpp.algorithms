
TEST(PairArrayLocal, simple) {
  using namespace imajuscule;

  LocalPairArray<int,float, 10> pa(2,4.f);

  EXPECT_EQ(10, pa.size());

  for(int i=0; i<10; ++i) {
    EXPECT_EQ(2, pa.firsts()[i]);
    EXPECT_FLOAT_EQ(4.f, pa.seconds()[i]);
  }
}


TEST(PairArrayLocal, AsBs) {
  using namespace imajuscule;
  LocalPairArray<BigAlignment,SmallAlignment, 10> pa(BigAlignment{1.0}, SmallAlignment{2.0});
  static_assert(decltype(pa)::order == detail::Order::As_Bs);


  EXPECT_EQ(10, pa.size());
  
  for(int i=0; i<10; ++i) {
    EXPECT_FLOAT_EQ(1.0, pa.firsts()[i].v);
    EXPECT_FLOAT_EQ(2.0, pa.seconds()[i].v);
  }
}

TEST(PairArrayLocal, BsAs) {
  using namespace imajuscule;
  LocalPairArray<SmallAlignment,BigAlignment, 10> pa(SmallAlignment{1.0}, BigAlignment{2.0});
  static_assert(decltype(pa)::order == detail::Order::Bs_As);

  EXPECT_EQ(10, pa.size());
  
  for(int i=0; i<10; ++i) {
    EXPECT_FLOAT_EQ(1.0, pa.firsts()[i].v);
    EXPECT_FLOAT_EQ(2.0, pa.seconds()[i].v);
  }
}

TEST(PairArrayLocal, RangeIteration) {
  using namespace imajuscule;
  LocalPairArray<int8_t,std::string, 10> pa(1, "hello");
  
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

TEST(PairArrayLocal, RangeIteration2) {
  using namespace imajuscule;
  LocalPairArray<std::string,int8_t, 10> pa("hello", 1);
  
  for(auto & p : firsts(pa)) {
    EXPECT_EQ("hello", p);
  }
  for(auto & p : seconds(pa)) {
    EXPECT_EQ(1, p);
  }
}
