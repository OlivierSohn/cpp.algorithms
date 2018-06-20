
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
