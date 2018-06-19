
TEST(PairArray, simple) {
  using namespace imajuscule;

  PairArray<int,float> pa(10,2,4.f);

  EXPECT_EQ(10, pa.size());

  for(int i=0; i<10; ++i) {
    EXPECT_EQ(2, pa.firsts()[i]);
    EXPECT_FLOAT_EQ(4.f, pa.seconds()[i]);
  }
}

namespace imajuscule {

  struct BigAlignment {
    BigAlignment() = default;
    BigAlignment(double d) {
      v = d;
    }
    union {
      std::aligned_storage_t<64, 64> placeholder;
      double v;
    };
  };

  struct SmallAlignment {
    double v;
  };
  
  static_assert(alignof(SmallAlignment)<alignof(BigAlignment));
  static_assert(alignof(SmallAlignment)>=alignof(void*));

}

TEST(PairArray, AsBs) {
  using namespace imajuscule;
  PairArray<BigAlignment,SmallAlignment> pa(10, BigAlignment{1.0}, SmallAlignment{2.0});
  static_assert(decltype(pa)::order == Order::As_Bs);


  EXPECT_EQ(10, pa.size());
  
  for(int i=0; i<10; ++i) {
    EXPECT_FLOAT_EQ(1.0, pa.firsts()[i].v);
    EXPECT_FLOAT_EQ(2.0, pa.seconds()[i].v);
  }
}

TEST(PairArray, BsAs) {
  using namespace imajuscule;
  PairArray<SmallAlignment,BigAlignment> pa(10, SmallAlignment{1.0}, BigAlignment{2.0});
  static_assert(decltype(pa)::order == Order::Bs_As);

  EXPECT_EQ(10, pa.size());
  
  for(int i=0; i<10; ++i) {
    EXPECT_FLOAT_EQ(1.0, pa.firsts()[i].v);
    EXPECT_FLOAT_EQ(2.0, pa.seconds()[i].v);
  }
}
