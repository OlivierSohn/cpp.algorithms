

TEST(Fifo1, cancelEmplace) {
  using namespace imajuscule;
  
  fifo1<int> f;

  f.emplace(3);
  EXPECT_FALSE(f.empty());
  EXPECT_TRUE(f.full());
  f.cancel_emplace();
  EXPECT_TRUE(f.empty());
  EXPECT_FALSE(f.full());
  
  f.emplace(3);
  EXPECT_FALSE(f.empty());
  EXPECT_TRUE(f.full());
  f.cancel_emplace();
  EXPECT_TRUE(f.empty());
  EXPECT_FALSE(f.full());
  
}

TEST(Fifo1, test) {
    using namespace imajuscule;
    fifo1<int> f;
    EXPECT_TRUE(f.empty());
    EXPECT_EQ(0,f.size());

    f.push(3);
    EXPECT_EQ(1,f.size());
    EXPECT_FALSE(f.empty());
    EXPECT_EQ(3,f.front());

    f.front() = 4;
    EXPECT_EQ(4,f.front());
    EXPECT_EQ(4,f.back());

    f.back()++;
    EXPECT_EQ(5,f.back());

    f.pop();
    EXPECT_TRUE(f.empty());
    EXPECT_EQ(0,f.size());

    int p = 6;
    f.push(p);
    EXPECT_EQ(1,f.size());
    EXPECT_FALSE(f.empty());
    EXPECT_EQ(6,f.front());

    f.pop();
    EXPECT_TRUE(f.empty());
    EXPECT_EQ(0,f.size());

    f.push(std::move(p));
    EXPECT_EQ(1,f.size());
    EXPECT_FALSE(f.empty());
    EXPECT_EQ(6,f.front());

    f.pop();
    EXPECT_TRUE(f.empty());
    EXPECT_EQ(0,f.size());

    f.emplace(7);
    EXPECT_EQ(1,f.size());
    EXPECT_FALSE(f.empty());
    EXPECT_EQ(7,f.front());

    fifo1<int> g;

    f.swap(g);
    EXPECT_TRUE(f.empty());
    EXPECT_EQ(0,f.size());

    EXPECT_EQ(1,g.size());
    EXPECT_FALSE(g.empty());
    EXPECT_EQ(7,g.front());
}
