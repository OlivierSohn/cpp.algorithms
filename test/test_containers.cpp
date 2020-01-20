

TEST(SlidingAverage, test) {
  using namespace imajuscule;
  slidingAverage<float> avg(3);
  EXPECT_FLOAT_EQ(0.f, avg.compute());
  EXPECT_FLOAT_EQ(0.f, avg.compute());
  avg.feed(3.f);
  EXPECT_EQ(3, avg.size());
  EXPECT_FLOAT_EQ(3.f, avg.compute());
  EXPECT_FLOAT_EQ(3.f, avg.compute());
  avg.feed(3.f);
  EXPECT_EQ(3, avg.size());
  EXPECT_FLOAT_EQ(3.f, avg.compute());
  EXPECT_FLOAT_EQ(3.f, avg.compute());
  avg.feed(0.f);
  EXPECT_EQ(3, avg.size());
  EXPECT_FLOAT_EQ(2.f, avg.compute());
  EXPECT_FLOAT_EQ(2.f, avg.compute());
  avg.feed(0.f);
  EXPECT_EQ(3, avg.size());
  EXPECT_FLOAT_EQ(1.f, avg.compute());
  EXPECT_FLOAT_EQ(1.f, avg.compute());
  avg.feed(0.f);
  EXPECT_EQ(3, avg.size());
  EXPECT_FLOAT_EQ(0.f, avg.compute());
  EXPECT_FLOAT_EQ(0.f, avg.compute());
  avg.feed(0.f);
  EXPECT_EQ(3, avg.size());
  EXPECT_FLOAT_EQ(0.f, avg.compute());
  EXPECT_FLOAT_EQ(0.f, avg.compute());
}


TEST(SlidingWindowAverage, test) {
    using namespace imajuscule;
    {
        slidingWindowedAverage<float> a(3, [](auto ratio){ return ratio; });
        a.feed(1);
        a.feed(1);
        a.feed(1);
        ASSERT_NEAR(2.f/3.f, a.compute(), 1e-6);
    }
    {
        slidingWindowedAverage<float> a(3, [](auto ratio){ return ratio; });
        a.feed(0);
        a.feed(1);
        a.feed(0);
        ASSERT_NEAR(1.f/3.f, a.compute(), 1e-6);
    }
    {
        slidingWindowedAverage<float> a(3, [](auto ratio){ return ratio; });
        a.feed(1);
        a.feed(0);
        a.feed(0);
        ASSERT_NEAR(.5f/3.f, a.compute(), 1e-6);
    }
    {
        slidingWindowedAverage<float> a(4, [](auto ratio){ return ratio; });
        a.feed(1);
        a.feed(1);
        a.feed(1);
        a.feed(1);
        ASSERT_NEAR(.6f, a.compute(), 1e-6);
    }
    {
        slidingWindowedAverage<float> a(4, [](auto ratio){ return ratio; });
        a.feed(0);
        a.feed(1);
        a.feed(1);
        a.feed(0);
        ASSERT_NEAR(.4f, a.compute(), 1e-6);
    }
    {
        slidingWindowedAverage<float> a(4, [](auto ratio){ return ratio; });
        a.feed(1);
        a.feed(0);
        a.feed(0);
        a.feed(1);
        ASSERT_NEAR(.2f, a.compute(), 1e-6);
    }
    {
        slidingWindowedAverage<float> a(4, [](auto ratio){ return ratio; });
        a.feed(0);
        a.feed(0);
        a.feed(0);
        a.feed(1);
        ASSERT_NEAR(.2f/2, a.compute(), 1e-6);
    }
    {
        slidingWindowedAverage<float> a(4, [](auto ratio){ return ratio; });
        a.feed(1);
        a.feed(0);
        a.feed(0);
        a.feed(0);
        ASSERT_NEAR(.2f/2, a.compute(), 1e-6);
    }
}

TEST(Container, ExtractFromEnd) {
  using namespace imajuscule;
  
  std::vector<std::unique_ptr<int>> v;
  auto p = std::make_unique<int>(5);
  auto ptr = p.get();
  v.push_back(std::move(p));
  int six;
  int * ptr2 = &six;
  EXPECT_EQ(nullptr, extractFromEnd(v,static_cast<int*>(0)).get());
  EXPECT_EQ(nullptr, extractFromEnd(v,ptr2).get());
  EXPECT_EQ(ptr    , extractFromEnd(v,ptr).get());
  EXPECT_TRUE(v.empty());
  EXPECT_EQ(nullptr, extractFromEnd(v,ptr).get());
}

TEST(Container, ExtractFromEnd_lasts) {
  using namespace imajuscule;

  std::vector<std::unique_ptr<int>> v;
  v.emplace_back(std::make_unique<int>(0));
  v.emplace_back(std::make_unique<int>(1));
  v.emplace_back(std::make_unique<int>(2));
  v.emplace_back(std::make_unique<int>(3));
  v.emplace_back(std::make_unique<int>(4));
  v.emplace_back(std::make_unique<int>(5));

  for(int i=5; i>=0; --i) {
    auto ptr = v[i].get();
    EXPECT_EQ(ptr, extractFromEnd(v,ptr).get());
  }
  EXPECT_TRUE(v.empty());
}

TEST(Container, ExtractFromEnd3_firsts) {
  using namespace imajuscule;
  
  std::vector<std::unique_ptr<int>> v;
  v.emplace_back(std::make_unique<int>(0));
  v.emplace_back(std::make_unique<int>(1));
  v.emplace_back(std::make_unique<int>(2));
  v.emplace_back(std::make_unique<int>(3));
  v.emplace_back(std::make_unique<int>(4));
  v.emplace_back(std::make_unique<int>(5));
  
  for(int i=5; i>=0; --i) {
    auto ptr = v[0].get();
    EXPECT_EQ(ptr, extractFromEnd(v,ptr).get());
  }
  EXPECT_TRUE(v.empty());
}

TEST(Container, ExtractFromEnd3_middle) {
  using namespace imajuscule;
  
  std::vector<std::unique_ptr<int>> v;
  v.emplace_back(std::make_unique<int>(0));
  v.emplace_back(std::make_unique<int>(1));
  v.emplace_back(std::make_unique<int>(2));
  v.emplace_back(std::make_unique<int>(3));
  v.emplace_back(std::make_unique<int>(4));
  v.emplace_back(std::make_unique<int>(5));
  
  for(int i=5; i>=0; --i) {
    auto ptr = v[i/2].get();
    EXPECT_EQ(ptr, extractFromEnd(v,ptr).get());
  }
  EXPECT_TRUE(v.empty());
}

struct TestLock {
  void lock() {
    locked = true;
  }
  void unlock() {
    locked = false;
  }
  bool locked = false;
};

TEST(Container, grow) {
  using namespace imajuscule;
  std::vector<int> v;
  auto w = mkVectorWrapper(v);
  TestLock l;
  
  EXPECT_TRUE(v.capacity() == 0);
  
  for(int i=0; i<10; ++i) {
    EXPECT_TRUE(reserveAndLock<CanRealloc::Yes>(0,w,l));
    EXPECT_TRUE(v.capacity() == 0);
    EXPECT_TRUE(l.locked);
    l.unlock();
  }
  
  for(int i=0; i<10; ++i) {
    EXPECT_TRUE(reserveAndLock<CanRealloc::Yes>(i,w,l));
    auto cap_now = v.capacity();
    EXPECT_LE(i,cap_now);
    EXPECT_TRUE(l.locked);
    l.unlock();
    
    EXPECT_TRUE(reserveAndLock<CanRealloc::Yes>(i,w,l));
    auto cap2 = v.capacity();
    EXPECT_EQ(cap2,cap_now);
    EXPECT_TRUE(l.locked);
    l.unlock();
    
  }
  
}

TEST(Container, grow2) {
  using namespace imajuscule;
  
  TestLock l;
  
  for(auto cap = 0; cap < 50; ++cap) {
    
    std::vector<int> v;
    auto w = mkVectorWrapper(v);

    v.resize(cap);
    EXPECT_TRUE(v.capacity() == cap);
    
    for(int i=0; i<100; ++i) {
      EXPECT_TRUE(reserveAndLock<CanRealloc::Yes>(i,w,l));
      auto cap_now = v.capacity();
      EXPECT_LE(cap+i,cap_now);
      EXPECT_TRUE(l.locked);
      l.unlock();
      
      EXPECT_TRUE(reserveAndLock<CanRealloc::Yes>(i,w,l));
      auto cap2 = v.capacity();
      EXPECT_EQ(cap2,cap_now);
      EXPECT_TRUE(l.locked);
      l.unlock();
    }
  }
  
}

TEST(Containers, fadeIn) {
  using namespace imajuscule;
  std::vector<float> c1 {12.f,12.f,12.f,12.f,12.f};
  
  auto c = withLinearFadeIn(3, c1);
  ASSERT_FLOAT_EQ(3.f, c[0]);
  ASSERT_FLOAT_EQ(6.f, c[1]);
  ASSERT_FLOAT_EQ(9.f, c[2]);
  ASSERT_FLOAT_EQ(12.f, c[3]);
  ASSERT_FLOAT_EQ(12.f, c[4]);
}

TEST(Containers, fadeOut) {
  using namespace imajuscule;
  std::vector<float> c1 {12.f,12.f,12.f,12.f,12.f};
  
  auto c = withLinearFadeOut(3, c1);
  ASSERT_FLOAT_EQ(12.f, c[0]);
  ASSERT_FLOAT_EQ(12.f, c[1]);
  ASSERT_FLOAT_EQ(9.f, c[2]);
  ASSERT_FLOAT_EQ(6.f, c[3]);
  ASSERT_FLOAT_EQ(3.f, c[4]);
}
