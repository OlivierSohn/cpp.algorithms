
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
