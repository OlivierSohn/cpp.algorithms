TEST(ForwardListLockfree, simple) {
  using namespace imajuscule::lockfree::scmp;
  forward_list<int> l;
  
  EXPECT_FALSE(l.try_pop_front());

  int n;
  auto inc = [&n] (auto &) {++n;};
  
  std::vector<int> v;
  auto record = [&v] (auto i) { v.push_back(i); };
  
  {
    n = 0;
    l.forEach(inc);
    EXPECT_EQ(0,n);
  }
  
  EXPECT_EQ(3,l.emplace_front(3));
  EXPECT_EQ(4,l.emplace_front(4));
  EXPECT_EQ(5,l.emplace_front(5));
  
  {
    n = 0;
    l.forEach(inc);
    EXPECT_EQ(3,n);
  }
  {
    n = 0;
    l.forEach(inc);
    EXPECT_EQ(3,n);
  }
  
  EXPECT_TRUE(l.try_pop_front());
  EXPECT_FALSE(l.try_pop_front());

  {
    n = 0;
    l.forEach(inc);
    EXPECT_EQ(2,n);
  }
  {
    n = 0;
    l.forEach(inc);
    EXPECT_EQ(2,n);
  }
  
  EXPECT_TRUE(l.try_pop_front());

  {
    n = 0;
    l.forEach(inc);
    EXPECT_EQ(1,n);
  }
  
  EXPECT_TRUE(l.try_pop_front());

  {
    n = 0;
    l.forEach(inc);
    EXPECT_EQ(0,n);
  }
  
  EXPECT_FALSE(l.try_pop_front());
  EXPECT_FALSE(l.try_pop_front());
  EXPECT_FALSE(l.try_pop_front());
  EXPECT_FALSE(l.try_pop_front());
  
  {
    n = 0;
    l.forEach(inc);
    EXPECT_EQ(0,n);
  }
  
  auto & r = l.emplace_front(8);
  
  {
    n = 0;
    l.forEach(inc);
    EXPECT_EQ(1,n);
  }
  
  {
    v.clear();
    l.forEach(record);
    EXPECT_EQ(std::vector<int>{8},v);
  }
  
  r = 7;
  
  {
    v.clear();
    l.forEach(record);
    EXPECT_EQ(std::vector<int>{7},v);
  }
  
  l.emplace_front(10);
  
  {
    v.clear();
    l.forEach(record);
    auto v2 = std::vector<int>{10,7};
    EXPECT_EQ(v2, v);
  }
  
}

TEST(ForwardListLockfree, garbage_pop1) {
  using namespace imajuscule::lockfree::scmp;
  using namespace imajuscule;
  
  EXPECT_EQ(0, Destructible::countLiveObjects());
  {
    forward_list<Destructible> l;
    EXPECT_EQ(0, Destructible::countLiveObjects());
    
    l.emplace_front();
    
    EXPECT_EQ(1, Destructible::countLiveObjects());
    
    EXPECT_TRUE(l.try_pop_front());

    EXPECT_EQ(1, Destructible::countLiveObjects());
    
    l.emplace_front();
    l.emplace_front();
    
    EXPECT_EQ(3, Destructible::countLiveObjects());
    
    int n;
    auto inc = [&n] (auto &) {++n;};
    
    {
      n = 0;
      l.forEach(inc);
      EXPECT_EQ(2,n); // one was removed
    }
    
    EXPECT_EQ(3, Destructible::countLiveObjects());
    
    l.emplace_front();
    
    // The removed elements is destroyed during the emplace_front call
    // following the forEach call following the removal.
    EXPECT_EQ(3, Destructible::countLiveObjects());
    
    l.emplace_front();
    l.emplace_front();
    
  }
  EXPECT_EQ(0, Destructible::countLiveObjects());
}


TEST(ForwardListLockfree, garbage_pop1_) {
  using namespace imajuscule::lockfree::scmp;
  using namespace imajuscule;
  
  EXPECT_EQ(0, Destructible::countLiveObjects());
  {
    forward_list<Destructible> l;
    EXPECT_EQ(0, Destructible::countLiveObjects());
    
    l.emplace_front();    
    EXPECT_EQ(1, Destructible::countLiveObjects());
    
    l.emplace_front();
    EXPECT_EQ(2, Destructible::countLiveObjects());

    l.emplace_front();
    EXPECT_EQ(3, Destructible::countLiveObjects());

    EXPECT_TRUE(l.try_pop_front());
    EXPECT_EQ(3, Destructible::countLiveObjects());
    
    int n;
    auto inc = [&n] (auto &) {++n;};
    
    {
      n = 0;
      l.forEach(inc);
      EXPECT_EQ(2,n); // one was removed
    }
    
    EXPECT_EQ(3, Destructible::countLiveObjects());
    
    l.emplace_front();
    
    // The removed elements is destroyed during the emplace_front call
    // following the forEach call following the removal.
    EXPECT_EQ(3, Destructible::countLiveObjects());
    
    l.emplace_front();
    l.emplace_front();
    
  }
  EXPECT_EQ(0, Destructible::countLiveObjects());
}


TEST(ForwardListLockfree, garbage_pop3) {
  using namespace imajuscule::lockfree::scmp;
  using namespace imajuscule;
  
  EXPECT_EQ(0, Destructible::countLiveObjects());
  {
    forward_list<Destructible> l;
    EXPECT_EQ(0, Destructible::countLiveObjects());
    
    l.emplace_front();
    l.emplace_front();
    l.emplace_front();
    
    EXPECT_EQ(3, Destructible::countLiveObjects());
    
    EXPECT_TRUE(l.try_pop_front());

    EXPECT_EQ(3, Destructible::countLiveObjects());
    
    int n;
    auto inc = [&n] (auto &) {++n;};
    
    {
      n = 0;
      l.forEach(inc);
      EXPECT_EQ(2,n); // one was removed
    }
    
    EXPECT_EQ(3, Destructible::countLiveObjects());
    
    l.emplace_front();
    
    // The removed element is destroyed during the emplace_front call
    // following the forEach call following the removal.
    EXPECT_EQ(3, Destructible::countLiveObjects());
    
    l.emplace_front();
    l.emplace_front();
    
  }
  EXPECT_EQ(0, Destructible::countLiveObjects());
}


TEST(ForwardListLockfree, garbage_popMid) {
  using namespace imajuscule::lockfree::scmp;
  using namespace imajuscule;
  
  EXPECT_EQ(0, Destructible::countLiveObjects());
  {
    forward_list<Destructible> l;
    EXPECT_EQ(0, Destructible::countLiveObjects());
    
    l.emplace_front();
    l.emplace_front();
    
    EXPECT_EQ(2, Destructible::countLiveObjects());
    
    EXPECT_TRUE(l.try_pop_front());

    EXPECT_EQ(2, Destructible::countLiveObjects());
    
    l.emplace_front();

    EXPECT_EQ(3, Destructible::countLiveObjects());

    int n;
    auto inc = [&n] (auto &) {++n;};
    
    {
      n = 0;
      l.forEach(inc);
      EXPECT_EQ(2,n); // one was removed
    }
    
    EXPECT_EQ(3, Destructible::countLiveObjects());
    
    l.emplace_front();
    
    // The removed element is destroyed during the emplace_front call
    // following the forEach call following the removal.
    EXPECT_EQ(3, Destructible::countLiveObjects());
    
    l.emplace_front();
    l.emplace_front();
    
  }
  EXPECT_EQ(0, Destructible::countLiveObjects());
}


TEST(ForwardListLockfree, garbage_popAntiMid) {
  using namespace imajuscule::lockfree::scmp;
  using namespace imajuscule;
  
  EXPECT_EQ(0, Destructible::countLiveObjects());
  {
    forward_list<Destructible> l;
    EXPECT_EQ(0, Destructible::countLiveObjects());
    
    l.emplace_front();
    EXPECT_EQ(1, Destructible::countLiveObjects());

    EXPECT_TRUE(l.try_pop_front());
    EXPECT_EQ(1, Destructible::countLiveObjects());

    l.emplace_front();
    EXPECT_EQ(2, Destructible::countLiveObjects());

    l.emplace_front();
    EXPECT_EQ(3, Destructible::countLiveObjects());

    EXPECT_TRUE(l.try_pop_front());
    EXPECT_EQ(3, Destructible::countLiveObjects());
    
    int n;
    auto inc = [&n] (auto &) {++n;};
    
    {
      n = 0;
      l.forEach(inc);
      EXPECT_EQ(1,n); // two were removed
    }
    
    EXPECT_EQ(3, Destructible::countLiveObjects());
    
    l.emplace_front();
    
    // The removed elements were destroyed during the emplace_front call
    // following the forEach call following the removal.
    EXPECT_EQ(2, Destructible::countLiveObjects());
    
    l.emplace_front();
    l.emplace_front();
    
  }
  EXPECT_EQ(0, Destructible::countLiveObjects());
}


TEST(ForwardListLockfree, garbage_pop4) {
  using namespace imajuscule::lockfree::scmp;
  using namespace imajuscule;
  
  EXPECT_EQ(0, Destructible::countLiveObjects());
  {
    forward_list<Destructible> l;
    EXPECT_EQ(0, Destructible::countLiveObjects());
    
    l.emplace_front();
    EXPECT_EQ(1, Destructible::countLiveObjects());
    
    EXPECT_TRUE(l.try_pop_front());
    EXPECT_EQ(1, Destructible::countLiveObjects());
    
    l.emplace_front();
    EXPECT_EQ(2, Destructible::countLiveObjects());
    
    EXPECT_TRUE(l.try_pop_front());
    EXPECT_EQ(2, Destructible::countLiveObjects());
    
    l.emplace_front();
    EXPECT_EQ(3, Destructible::countLiveObjects());
    
    EXPECT_TRUE(l.try_pop_front());
    EXPECT_EQ(3, Destructible::countLiveObjects());
    
    l.emplace_front();
    EXPECT_EQ(4, Destructible::countLiveObjects());
    
    EXPECT_TRUE(l.try_pop_front());
    EXPECT_EQ(4, Destructible::countLiveObjects());

    int n;
    auto inc = [&n] (auto &) {++n;};
    
    {
      n = 0;
      l.forEach(inc);
      EXPECT_EQ(0,n);
    }
    
    EXPECT_EQ(4, Destructible::countLiveObjects());
    
    l.emplace_front();
    
    // The removed elements were destroyed during the emplace_front call
    // following the forEach call following the removal.
    EXPECT_EQ(1, Destructible::countLiveObjects());
  }
  EXPECT_EQ(0, Destructible::countLiveObjects());
}


TEST(ForwardListLockfree, garbage_popAll) {
  using namespace imajuscule::lockfree::scmp;
  using namespace imajuscule;
  
  EXPECT_EQ(0, Destructible::countLiveObjects());
  {
    forward_list<Destructible> l;
    EXPECT_EQ(0, Destructible::countLiveObjects());
    
    l.emplace_front();
    EXPECT_EQ(1, Destructible::countLiveObjects());
    
    EXPECT_TRUE(l.try_pop_front());
    EXPECT_EQ(1, Destructible::countLiveObjects());
    
    l.emplace_front();
    EXPECT_EQ(2, Destructible::countLiveObjects());

    EXPECT_TRUE(l.try_pop_front());
    EXPECT_EQ(2, Destructible::countLiveObjects());

    l.emplace_front();
    EXPECT_EQ(3, Destructible::countLiveObjects());
    
    EXPECT_TRUE(l.try_pop_front());
    EXPECT_EQ(3, Destructible::countLiveObjects());
    
    int n;
    auto inc = [&n] (auto &) {++n;};
    
    {
      n = 0;
      l.forEach(inc);
      EXPECT_EQ(0,n); // 3 were removed
    }
    
    EXPECT_EQ(3, Destructible::countLiveObjects());
    
    l.emplace_front();
    
    // The removed elements were destroyed during the emplace_front call
    // following the forEach call following the removal.
    EXPECT_EQ(1, Destructible::countLiveObjects());
    
    l.emplace_front();
    l.emplace_front();
    
  }
  EXPECT_EQ(0, Destructible::countLiveObjects());
}


TEST(ForwardListLockfree, consumerProducer) {
  using namespace imajuscule::lockfree::scmp;

  std::atomic_bool run = true;
  forward_list<int> l;

  bool ok = true;
  std::thread produce([&ok, &l, &run](){
    
    int n=0;
    while(run) {
      l.emplace_front(n);
      // this 'try_pop_front' is guaranteed to succeed because it is the first
      // 'try_pop_front' since the last 'emplace_front'
      ok = l.try_pop_front() && ok;
      ++n;
    }
  });
  
  std::vector<int> read;
  read.reserve(1000000);

  std::thread schedule([&run](){
    std::this_thread::sleep_for(std::chrono::milliseconds(1000));
    run = false;
  });

  
  while(run) {
    l.forEach([&read](int i) { read.push_back(i); });
  }

  produce.join();
  schedule.join();

  EXPECT_TRUE(ok);
}


TEST(ForwardListLockfree, consumerProducerBig) {
  using namespace imajuscule::lockfree::scmp;
  using namespace imajuscule;
  
  std::atomic_bool run = true;
  forward_list<RepeatInt> l;
  
  bool ok = true;
  std::thread produce([&ok, &l, &run](){
    int n=0;
    while(run) {
      l.emplace_front(n);
      // this 'try_pop_front' is guaranteed to succeed because it is the first
      // 'try_pop_front' since the last 'emplace_front'
      ok = l.try_pop_front() && ok;
      ++n;
    }
  });
  
  std::vector<int> read;
  read.reserve(1000000);
  
  std::thread schedule([&run](){
    std::this_thread::sleep_for(std::chrono::milliseconds(1000));
    run = false;
  });
  
  
  while(run) {
    l.forEach([&read](auto & i) {
      EXPECT_TRUE(i.isValid());
      read.push_back(i.get());
    });
  }
  
  produce.join();
  schedule.join();
  
  EXPECT_TRUE(ok);
}


TEST(ForwardListLockfree, consumerProducers) {
  using namespace imajuscule::lockfree::scmp;
  
  std::atomic_bool run = true;
  forward_list<int> l;
  
  constexpr auto nthreads = 20;
  
  std::vector<std::thread> v_threads;
  v_threads.reserve(nthreads);
  
  for(int i=0; i<nthreads; ++i) {
    // producers
    v_threads.emplace_back(std::thread{[&l, &run](){
      int n=0;
      while(run) {
        l.emplace_front(n);
        // this 'try_pop_front' is not guaranteed to succeed
        // because it may not be the first
        // 'try_pop_front' since the last 'emplace_front'
        l.try_pop_front();
        ++n;
      }
    }});
  }

  std::vector<int> read;
  read.reserve(1000000);
  
  std::thread schedule([&run](){
    std::this_thread::sleep_for(std::chrono::milliseconds(1000));
    run = false;
  });
  
  
  while(run) {
    l.forEach([&read](int i) { read.push_back(i); });
  }
  
  for(auto & t: v_threads) {
    t.join();
  }
  schedule.join();
  
}


TEST(ForwardListLockfree, consumerProducersBig) {
  using namespace imajuscule::lockfree::scmp;
  using namespace imajuscule;
  
  std::atomic_bool run = true;
  forward_list<RepeatInt> l;
  
  constexpr auto nthreads = 20;
  
  std::vector<std::thread> v_threads;
  v_threads.reserve(nthreads);
  
  for(int i=0; i<nthreads; ++i) {
    // producers
    v_threads.emplace_back(std::thread{[&l, &run](){
      int n=0;
      while(run) {
        l.emplace_front(n);
        // this 'try_pop_front' is not guaranteed to succeed
        // because it may not be the first
        // 'try_pop_front' since the last 'emplace_front'
        l.try_pop_front();
        ++n;
      }
    }});
  }
  
  std::vector<int> read;
  read.reserve(1000000);
  
  std::thread schedule([&run](){
    std::this_thread::sleep_for(std::chrono::milliseconds(1000));
    run = false;
  });
  
  
  while(run) {
    l.forEach([&read](auto & i) {
      EXPECT_TRUE(i.isValid());
      read.push_back(i.get());
    });
  }
  
  for(auto & t: v_threads) {
    t.join();
  }
  schedule.join();
  
}
