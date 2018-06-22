TEST(FifoLockFreeSCMP, singleThread_0_elt) {
  using namespace imajuscule;
  
  lockfree::scmp::fifo<int> a(0);

  EXPECT_FALSE(a.tryDequeue([](auto &){}));
  EXPECT_FALSE(a.tryEnqueue(1));
  EXPECT_FALSE(a.tryDequeue([](auto &){}));
}

TEST(FifoLockFreeSCMP, singleThread_1_elt) {
  using namespace imajuscule;
  
  lockfree::scmp::fifo<int> a(1);
  
  EXPECT_FALSE(a.tryDequeue([](auto &){}));
  EXPECT_TRUE(a.tryEnqueue(1));
  {
    int v=0;
    auto d = a.tryDequeue([&v](int i){v = i;});
    EXPECT_TRUE(d);
    EXPECT_EQ(1,v);
  }

  EXPECT_TRUE(a.tryEnqueue(2));
  {
    int v=0;
    auto d = a.tryDequeue([&v](int i){v = i;});
    EXPECT_TRUE(d);
    EXPECT_EQ(2,v);
  }

  EXPECT_TRUE(a.tryEnqueue(3));
  EXPECT_FALSE(a.tryEnqueue(4));
  {
    int v=0;
    auto d = a.tryDequeue([&v](int i){v = i;});
    EXPECT_TRUE(d);
    EXPECT_EQ(3,v);
  }

  EXPECT_TRUE(a.tryEnqueue(5));
  {
    int v=0;
    auto d = a.tryDequeue([&v](int i){v = i;});
    EXPECT_TRUE(d);
    EXPECT_EQ(5,v);
  }
  EXPECT_FALSE(a.tryDequeue([](auto &){}));
}

TEST(FifoLockFreeSCMP, singleThread_2_elts) {
  using namespace imajuscule;
  
  lockfree::scmp::fifo<int> a(2);

  EXPECT_FALSE(a.tryDequeue([](auto &){}));
  EXPECT_TRUE(a.tryEnqueue(1));
  {
    int v=0;
    auto d = a.tryDequeue([&v](int i){v = i;});
    EXPECT_TRUE(d);
    EXPECT_EQ(1,v);
  }
  EXPECT_TRUE(a.tryEnqueue(2));
  {
    int v=0;
    auto d = a.tryDequeue([&v](int i){v = i;});
    EXPECT_TRUE(d);
    EXPECT_EQ(2,v);
  }
  EXPECT_TRUE(a.tryEnqueue(3));
  {
    int v=0;
    auto d = a.tryDequeue([&v](int i){v = i;});
    EXPECT_TRUE(d);
    EXPECT_EQ(3,v);
  }
  EXPECT_FALSE(a.tryDequeue([](auto &){}));

  EXPECT_TRUE(a.tryEnqueue(4));
  EXPECT_TRUE(a.tryEnqueue(5));
  EXPECT_FALSE(a.tryEnqueue(6));
  {
    int v=0;
    auto d = a.tryDequeue([&v](int i){v = i;});
    EXPECT_TRUE(d);
    EXPECT_EQ(4,v);
  }
  {
    int v=0;
    auto d = a.tryDequeue([&v](int i){v = i;});
    EXPECT_TRUE(d);
    EXPECT_EQ(5,v);
  }
  EXPECT_FALSE(a.tryDequeue([](auto &){}));

  EXPECT_TRUE(a.tryEnqueue(7));
  EXPECT_TRUE(a.tryEnqueue(8));
  EXPECT_FALSE(a.tryEnqueue(9));
  {
    int v=0;
    auto d = a.tryDequeue([&v](int i){v = i;});
    EXPECT_TRUE(d);
    EXPECT_EQ(7,v);
  }
  EXPECT_TRUE(a.tryEnqueue(10));
  {
    int v=0;
    auto d = a.tryDequeue([&v](int i){v = i;});
    EXPECT_TRUE(d);
    EXPECT_EQ(8,v);
  }
  {
    int v=0;
    auto d = a.tryDequeue([&v](int i){v = i;});
    EXPECT_TRUE(d);
    EXPECT_EQ(10,v);
  }
  EXPECT_FALSE(a.tryDequeue([](auto &){}));
}

template<int arraySz, typename T>
void testLFFIFOSlowConsumer() {
  using namespace imajuscule;
  
  LG(INFO,"using array size %d", arraySz);
  
  lockfree::scmp::fifo<T> a(arraySz);
  
  constexpr auto nthreads = 20;
  
  std::vector<std::thread> v_threads;
  v_threads.reserve(nthreads);
  std::atomic<int> nStartedProducers = 0;
  std::atomic<int> nStoppedProducers = 0;

#ifdef NDEBUG
  constexpr auto nNumbers = 10000;
#else
  constexpr auto nNumbers = 100;
#endif

  std::vector<int> values(nNumbers);
  std::iota(values.begin(), values.end(), 0);
  
  for(int i=0; i<nthreads; ++i) {
    // These are the producers : they add values to the array in a random order.
    v_threads.emplace_back(std::thread{[&a,&values,&nStartedProducers,&nStoppedProducers] () {
      auto v = values;
      Shuffle(v);
      ++nStartedProducers;
      while(nStartedProducers != nthreads);
      for(auto e : v) {
        while(!a.tryEnqueue({e})); // repeat until we succeed.
      }
      ++nStoppedProducers;
    }});
  }
  
  // I am the consumer
  std::vector<int> received;
  std::vector<bool> valid;
  
  received.reserve(nNumbers * nthreads);
  valid.reserve(nNumbers * nthreads);
  
  int nDequeued = 0;
  
  int nConsecutiveFails = 0;
  bool oneMoreChance = true;
  while(received.size() < received.capacity()) {
    std::this_thread::sleep_for(std::chrono::microseconds(10));

    auto nStoppedProd = nStoppedProducers.load();
    
    auto n = a.dequeueAll([&valid, &received](auto & ri){
      if constexpr (std::is_same_v<T, RepeatInt>) {
        valid.push_back(ri.isValid());
        received.push_back(ri.get());
      }
      else {
        valid.push_back(true);
        received.push_back(ri);
      }
    });
    
    nDequeued += n;
    if(n) {
      nConsecutiveFails = 0;
    }
    else {
      ++nConsecutiveFails;
      if(nStoppedProd == nthreads) {
        if(oneMoreChance) {
          oneMoreChance = false;
          continue;
        }
        LG(WARN,"aborting test, %d consecutive fails with %d stopped producers. %d / %d",
           nConsecutiveFails,
           nStoppedProd,
           received.size(),
           received.capacity());
        EXPECT_FALSE(true);
        break;
      }
    }
    if(nConsecutiveFails == 1000000) {
      LG(WARN,"%d consecutive fails with %d stopped producers. %d / %d",
         nConsecutiveFails,
         nStoppedProd,
         received.size(),
         received.capacity());
    }
  }

  EXPECT_EQ(nDequeued, received.size());
  
  for(auto & t: v_threads) {
    t.join();
  }
  
  for(auto b : valid) {
    ASSERT_TRUE(b);
  }
  
  std::vector<int> expected;
  expected.reserve(nNumbers * nthreads);
  for(int i=0; i<nthreads;++i) {
    expected.insert(expected.end(), values.begin(), values.end());
  }
  
  StdSort(expected);
  StdSort(received);
  EXPECT_EQ(expected,received);
}


template<int arraySz, typename T>
void testLFFIFOFastConsumer() {
  using namespace imajuscule;
  
  LG(INFO,"using array size %d and scheduling", arraySz);
  
  lockfree::scmp::fifo<T> a(arraySz);
  
  constexpr auto nthreads = 20;
  
  std::vector<std::thread> v_threads;
  v_threads.reserve(nthreads);
  
#ifdef NDEBUG
  constexpr auto nNumbers = 10000;
#else
  constexpr auto nNumbers = 100;
#endif
  std::vector<int> values(nNumbers);
  std::iota(values.begin(), values.end(), 0);
  std::atomic<int> nStartedProducers = 0;
  std::atomic<int> nStoppedProducers = 0;
  std::atomic<bool> run = false;
  std::atomic<bool> run2 = true;
  
  for(int i=0; i<nthreads; ++i) {
    // These are the producers : they add values to the array in a random order.
    v_threads.emplace_back(std::thread{[&a,&values,&run,&nStartedProducers,&nStoppedProducers] () {
      ++nStartedProducers;
      auto v = values;
      Shuffle(v);

      while(nStartedProducers != nthreads);
      for(auto e : v) {
        while(!(run && a.tryEnqueue({e}))); // repeat until we succeed.
      }
      ++nStoppedProducers;
    }});
  }
  
  std::thread sched([&run, &run2] () {
    while(run2) {
      run = true;
      std::this_thread::sleep_for(std::chrono::milliseconds(10));
      run = false;
      std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
  });
  
  // I am the consumer
  std::vector<int> received;
  std::vector<bool> valid;
  
  received.reserve(nNumbers * nthreads);
  valid.reserve(nNumbers * nthreads);
  
  int nDequeued = 0;

  int nConsecutiveFails = 0;
  bool oneMoreChance = true;
  while(received.size() < received.capacity()) {
    auto nStoppedProd = nStoppedProducers.load();

    auto n = a.dequeueAll([&valid, &received](auto & ri){
      if constexpr (std::is_same_v<T, RepeatInt>) {
        valid.push_back(ri.isValid());
        received.push_back(ri.get());
      }
      else {
        valid.push_back(true);
        received.push_back(ri);
      }
    });
    
    nDequeued += n;
    
    if(!run || n) {
      nConsecutiveFails = 0;
    }
    else {
      ++nConsecutiveFails;
      if(nStoppedProd == nthreads) {
        if(oneMoreChance) {
          oneMoreChance = false;
          continue;
        }
        LG(WARN,"aborting test, %d consecutive fails with %d stopped producers. %d / %d",
           nConsecutiveFails,
           nStoppedProd,
           received.size(),
           received.capacity());
        EXPECT_FALSE(true);
        break;
      }
    }
    if(nConsecutiveFails == 1000000) {
      LG(WARN,"%d consecutive fails with %d stopped producers. %d / %d",
         nConsecutiveFails,
         nStoppedProd,
         received.size(),
         received.capacity());
    }
  }
  
  EXPECT_EQ(nDequeued, received.size());
  
  for(auto & t: v_threads) {
    t.join();
  }
  run2 = false;
  sched.join();
  
  for(auto b : valid) {
    ASSERT_TRUE(b);
  }
  
  std::vector<int> expected;
  expected.reserve(nNumbers * nthreads);
  for(int i=0; i<nthreads;++i) {
    expected.insert(expected.end(), values.begin(), values.end());
  }
  
  StdSort(expected);
  StdSort(received);
  EXPECT_EQ(expected,received);
}

template<int arraySz>
void testLFFIFO() {
  using namespace imajuscule;
  testLFFIFOFastConsumer<arraySz, int>();
  testLFFIFOFastConsumer<arraySz, RepeatInt>();

  testLFFIFOSlowConsumer<arraySz, int>();
  testLFFIFOSlowConsumer<arraySz, RepeatInt>();
}

TEST(FifoLockFreeSCMP, multiThread) {
  testLFFIFO<1>();
  testLFFIFO<2>();
  testLFFIFO<3>();
  testLFFIFO<10>();
  testLFFIFO<100>();
  testLFFIFO<1000>();
  testLFFIFO<10000>();
  testLFFIFO<100000>();
}

TEST(FifoLockFreeSCMP, OnRemovalAssignFromDefault_UP) {
  using namespace imajuscule;
  
  EXPECT_EQ(0,nLiveObjects());
  {
    lockfree::scmp::fifo<std::unique_ptr<Destructible>, lockfree::scmp::OnRemoval::AssignFromDefault> a(10);
    
    a.tryEnqueue(std::make_unique<Destructible>());
    EXPECT_EQ(1,nLiveObjects());
    a.tryEnqueue(std::make_unique<Destructible>());
    EXPECT_EQ(2,nLiveObjects());
    
    // remove all elements
    a.dequeueAll([](auto &){});

    EXPECT_EQ(0,nLiveObjects());
  }
  EXPECT_EQ(0,nLiveObjects());
}

TEST(FifoLockFreeSCMP, OnRemovalAssignFromDefault) {
  using namespace imajuscule;

  EXPECT_EQ(0,nLiveObjects());
  {
    lockfree::scmp::fifo<Destructible, lockfree::scmp::OnRemoval::AssignFromDefault> a(10);
    EXPECT_EQ(11,nLiveObjects());

    a.tryEnqueue({});
    EXPECT_EQ(11,nLiveObjects());
    a.tryEnqueue({});
    EXPECT_EQ(11,nLiveObjects());
    
    // remove all elements
    a.dequeueAll([](auto &){});

    EXPECT_EQ(11,nLiveObjects());
  }
  EXPECT_EQ(0,nLiveObjects());
}


TEST(FifoLockFreeSCMP, OnRemovalDoNothing_UP) {
  using namespace imajuscule;
  
  EXPECT_EQ(0,nLiveObjects());
  {
    lockfree::scmp::fifo<std::unique_ptr<Destructible>, lockfree::scmp::OnRemoval::DoNothing> a(10);
    
    a.tryEnqueue(std::make_unique<Destructible>());
    EXPECT_EQ(1,nLiveObjects());
    a.tryEnqueue(std::make_unique<Destructible>());
    EXPECT_EQ(2,nLiveObjects());
    
    // remove all elements
    a.dequeueAll([](auto &){});
    EXPECT_EQ(2,nLiveObjects()); // the removal hasn't changed the number of live objects

    // the next 9 inserts will succeed and change the number of live objects
    EXPECT_TRUE(a.tryEnqueue(std::make_unique<Destructible>()));
    EXPECT_EQ(3,nLiveObjects());
    EXPECT_TRUE(a.tryEnqueue(std::make_unique<Destructible>()));
    EXPECT_EQ(4,nLiveObjects());
    EXPECT_TRUE(a.tryEnqueue(std::make_unique<Destructible>()));
    EXPECT_EQ(5,nLiveObjects());
    EXPECT_TRUE(a.tryEnqueue(std::make_unique<Destructible>()));
    EXPECT_EQ(6,nLiveObjects());
    EXPECT_TRUE(a.tryEnqueue(std::make_unique<Destructible>()));
    EXPECT_EQ(7,nLiveObjects());
    EXPECT_TRUE(a.tryEnqueue(std::make_unique<Destructible>()));
    EXPECT_EQ(8,nLiveObjects());
    EXPECT_TRUE(a.tryEnqueue(std::make_unique<Destructible>()));
    EXPECT_EQ(9,nLiveObjects());
    EXPECT_TRUE(a.tryEnqueue(std::make_unique<Destructible>()));
    EXPECT_EQ(10,nLiveObjects());
    EXPECT_TRUE(a.tryEnqueue(std::make_unique<Destructible>()));
    EXPECT_EQ(11,nLiveObjects());

    // the next inserts will succeed and not change the number of live objects
    EXPECT_TRUE(a.tryEnqueue(std::make_unique<Destructible>()));
    EXPECT_EQ(11,nLiveObjects());
  }
  EXPECT_EQ(0,nLiveObjects());
}

TEST(FifoLockFreeSCMP, OnRemovalDoNothing) {
  using namespace imajuscule;
  
  EXPECT_EQ(0,nLiveObjects());
  {
    lockfree::scmp::fifo<Destructible, lockfree::scmp::OnRemoval::DoNothing> a(10);
    EXPECT_EQ(11,nLiveObjects());

    a.tryEnqueue({});
    EXPECT_EQ(11,nLiveObjects());
    a.tryEnqueue({});
    EXPECT_EQ(11,nLiveObjects());
    
    // remove all elements
    a.dequeueAll([](auto &){});
    EXPECT_EQ(11,nLiveObjects());

    a.tryEnqueue({});
    EXPECT_EQ(11,nLiveObjects());
    a.tryEnqueue({});
    EXPECT_EQ(11,nLiveObjects());
  }
  EXPECT_EQ(0,nLiveObjects());
}


