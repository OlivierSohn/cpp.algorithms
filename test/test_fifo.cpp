
TEST(Fifo, cancelEmplace) {
  using namespace imajuscule;
  
  fifo<int> f(3);
  f.emplace(3);
  EXPECT_FALSE(f.empty());
  f.cancel_emplace();
  EXPECT_TRUE(f.empty());

  f.emplace(3);
  f.emplace(3);
  f.emplace(3);
  EXPECT_TRUE(f.full());

  f.cancel_emplace();
  EXPECT_FALSE(f.empty());
  f.cancel_emplace();
  EXPECT_FALSE(f.empty());
  f.cancel_emplace();
  EXPECT_TRUE(f.empty());
}


TEST(Fifo, unsafeOneElt) {
    using namespace imajuscule;

    fifo<int> f(1);
    EXPECT_TRUE(f.empty());
    EXPECT_FALSE(f.full());
    EXPECT_EQ(0,f.size());

    f.emplace(3);
    EXPECT_EQ(1,f.size());
    EXPECT_FALSE(f.empty());
    EXPECT_TRUE(f.full());
    EXPECT_EQ(3,f.front());

    f.front() = 4;
    EXPECT_EQ(4,f.front());
    EXPECT_EQ(4,f.back());

    f.back()++;
    EXPECT_EQ(5,f.back());

    f.pop();
    EXPECT_TRUE(f.empty());
    EXPECT_FALSE(f.full());
    EXPECT_EQ(0,f.size());

    int p = 6;
    f.emplace(int{p});
    EXPECT_EQ(1,f.size());
    EXPECT_FALSE(f.empty());
    EXPECT_TRUE(f.full());
    EXPECT_EQ(6,f.front());

    f.pop();
    EXPECT_TRUE(f.empty());
    EXPECT_FALSE(f.full());
    EXPECT_EQ(0,f.size());

    f.emplace(int{p});
    EXPECT_EQ(1,f.size());
    EXPECT_FALSE(f.empty());
    EXPECT_TRUE(f.full());
    EXPECT_EQ(6,f.front());

    f.pop();
    EXPECT_TRUE(f.empty());
    EXPECT_FALSE(f.full());
    EXPECT_EQ(0,f.size());

    f.emplace(7);
    EXPECT_EQ(1,f.size());
    EXPECT_FALSE(f.empty());
    EXPECT_TRUE(f.full());
    EXPECT_EQ(7,f.front());

    fifo<int> g;

    f.swap(g);
    EXPECT_TRUE(f.empty());
    EXPECT_FALSE(f.full());
    EXPECT_EQ(0,f.size());

    EXPECT_EQ(1,g.size());
    EXPECT_FALSE(g.empty());
    EXPECT_TRUE(g.full());
    EXPECT_EQ(7,g.front());
}

TEST(Fifo, testMultipleEltSimple) {
  using namespace imajuscule;
  
  constexpr auto initialCapacity = 54;
  constexpr auto finalCapacity = 154;
  
  fifo<int> f(initialCapacity);
  
  EXPECT_TRUE(f.empty());
  EXPECT_FALSE(f.full());
  EXPECT_EQ(0,f.size());
  
  for(int i=0; i<initialCapacity; ++i) {
    f.emplace(int{i});
    
    EXPECT_EQ(i+1,f.size());
    EXPECT_FALSE(f.empty());
    EXPECT_EQ(0,f.front());
    EXPECT_EQ(i,f.back());
    EXPECT_EQ(i==initialCapacity-1, f.full());
  }

  // by now, the queue contains consecutive ints from 0 to (initialCapacity-1)

  std::vector<int> s;
  s.reserve(finalCapacity + 1); // + 1
  
  EXPECT_TRUE(f.full());
  
  f.trySwapUnderlyingContainer(s);
  
  EXPECT_FALSE(f.full());
  
  for(int i=initialCapacity; i<finalCapacity; ++i) {
    f.emplace(int{i});
    
    EXPECT_EQ(i+1,f.size());
    EXPECT_FALSE(f.empty());
    EXPECT_EQ(0,f.front());
    EXPECT_EQ(i,f.back());
    EXPECT_EQ(i==finalCapacity-1, f.full());
  }
  
  for(int i=0; i<finalCapacity; ++i) {
    EXPECT_EQ(i, f.front());
    EXPECT_FALSE(f.empty());
    
    f.pop();
    
    EXPECT_FALSE(f.full());
  }
  EXPECT_TRUE(f.empty());
}

TEST(Fifo, testMultipleEltConsumeInbetween) {
  using namespace imajuscule;
  
  constexpr auto initialCapacity = 54;
  constexpr auto finalCapacity = 154;

  constexpr auto consume = 10;

  fifo<int> f(initialCapacity);
  
  EXPECT_TRUE(f.empty());
  EXPECT_FALSE(f.full());
  EXPECT_EQ(0,f.size());
  
  for(int i=0; i<initialCapacity; ++i) {
    f.emplace(i-consume);
    
    EXPECT_EQ(i+1,f.size());
    EXPECT_FALSE(f.empty());
    EXPECT_EQ(0-consume,f.front());
    EXPECT_EQ(i-consume,f.back());
    EXPECT_EQ(i==initialCapacity-1, f.full());
  }

  for(int i=0; i<consume; ++i) {
    EXPECT_EQ(i-consume, f.front());
    EXPECT_FALSE(f.empty());
    
    f.pop();
    
    EXPECT_FALSE(f.full());
  }

  EXPECT_EQ(0,f.front());

  for(int i=0; i<consume; ++i) {
    f.emplace(i+initialCapacity-consume);
    
    EXPECT_EQ(i+initialCapacity-consume+1,f.size());
    EXPECT_FALSE(f.empty());
    EXPECT_EQ(0,f.front());
    EXPECT_EQ(i+initialCapacity-consume,f.back());
    EXPECT_EQ(i==consume-1, f.full());
  }
  
  // by now, the queue contains consecutive ints from 0 to (initialCapacity-1)

  std::vector<int> s;
  s.reserve(finalCapacity + 1); // + 1
  
  EXPECT_TRUE(f.full());
  
  f.trySwapUnderlyingContainer(s);
  
  EXPECT_FALSE(f.full());
  
  for(int i=initialCapacity; i<finalCapacity; ++i) {
    f.emplace(int{i});
    
    EXPECT_EQ(i+1,f.size());
    EXPECT_FALSE(f.empty());
    EXPECT_EQ(0,f.front());
    EXPECT_EQ(i,f.back());
    EXPECT_EQ(i==finalCapacity-1, f.full());
  }
  
  for(int i=0; i<finalCapacity; ++i) {
    EXPECT_EQ(i, f.front());
    EXPECT_FALSE(f.empty());
    
    f.pop();
    
    EXPECT_FALSE(f.full());
  }
  EXPECT_TRUE(f.empty());
}

TEST(Fifo, testNoRealloc) {
  using namespace imajuscule;
  
  fifo<int> f(0);
  std::mutex m;
  EXPECT_FALSE(safeTryEmplace<CanRealloc::No>(f, m, 48));
  EXPECT_TRUE(f.empty());
}

TEST(Fifo, testLock) {
  using namespace imajuscule;
  
  fifo<int> f;
  std::mutex m;
  
  constexpr auto nElems = 1000;
  for(int i=0; i<nElems; ++i) {
    EXPECT_TRUE(safeTryEmplace<CanRealloc::Yes>(f, m, std::move(i)));
  }
  
  for(int i=0; i<nElems; ++i) {
    EXPECT_EQ(i, f.front());
    EXPECT_FALSE(f.empty());
    
    f.pop();
    
    EXPECT_FALSE(f.full());
  }
  EXPECT_TRUE(f.empty());
}
