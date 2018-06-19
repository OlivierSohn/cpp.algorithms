TEST(StaticVectorLockFreeSCMP, singleThread_0_elt) {
  using namespace imajuscule;
  
  lockfree::scmp::static_vector<int> a(0);
  
  {
    int n = 0;
    int nRemoved = a.forEach([&n](int & elt) {
      ++n;
      return true;
    });
    EXPECT_EQ(0,nRemoved);
    EXPECT_EQ(0,n);
  }
  
  EXPECT_FALSE(a.tryInsert(1));
  
  {
    int n = 0;
    int nRemoved = a.forEach([&n](int & elt) {
      ++n;
      return true;
    });
    EXPECT_EQ(0,nRemoved);
    EXPECT_EQ(0,n);
  }
}

TEST(StaticVectorLockFreeSCMP, singleThread_1_elt) {
  using namespace imajuscule;
  
  lockfree::scmp::static_vector<int> a(1);
  
  {
    int n = 0;
    int nRemoved = a.forEach([&n](int & elt) {
      ++n;
      return true;
    });
    EXPECT_EQ(0,nRemoved);
    EXPECT_EQ(0,n);
  }
  
  EXPECT_TRUE(a.tryInsert(1)); // add an element
  
  {
    int n = 0;
    int nRemoved = a.forEach([&n](int & elt) {
      if(elt == 1) {
        elt = 8;
        ++n;
        return true;
      }
      return false;
    });
    EXPECT_EQ(0,nRemoved);
    EXPECT_EQ(1,n);
  }
  
  {
    int n = 0;
    int nRemoved = a.forEach([&n](int & elt) {
      if(elt == 8) {
        ++n;
        return false; // delete the element
      }
      return true;
    });
    EXPECT_EQ(1,nRemoved);
    EXPECT_EQ(1,n);
  }
  
  {
    int n = 0;
    int nRemoved = a.forEach([&n](int & elt) {
      ++n;
      return true;
    });
    EXPECT_EQ(0,nRemoved);
    EXPECT_EQ(0,n); // because the element was deleted
  }
  
  EXPECT_TRUE(a.tryInsert(4)); // add an element
  EXPECT_FALSE(a.tryInsert(5)); // add an element
  EXPECT_FALSE(a.tryInsert(6)); // add an element
  EXPECT_FALSE(a.tryInsert(7)); // add an element
  
  {
    int n = 0;
    int nRemoved = a.forEach([&n](int & elt) {
      if(elt == 4) {
        ++n;
        return true;
      }
      return false;
    });
    EXPECT_EQ(0,nRemoved);
    EXPECT_EQ(1,n);
  }
  
}

TEST(StaticVectorLockFreeSCMP, singleThread_2_elts) {
  using namespace imajuscule;
  
  lockfree::scmp::static_vector<int> a(2);
  
  {
    int n = 0;
    int nRemoved = a.forEach([&n](int & elt) {
      ++n;
      return true;
    });
    EXPECT_EQ(0,nRemoved);
    EXPECT_EQ(0,n);
  }
  
  EXPECT_TRUE(a.tryInsert(1)); // add an element
  
  {
    int n = 0;
    int nRemoved = a.forEach([&n](int & elt) {
      if(elt == 1) {
        elt = 8;
        ++n;
        return true;
      }
      return false;
    });
    EXPECT_EQ(0,nRemoved);
    EXPECT_EQ(1,n);
  }
  
  EXPECT_TRUE(a.tryInsert(2)); // add another element
  
  {
    int n = 0;
    int nRemoved = a.forEach([&n](int & elt) {
      if(elt == 2) {
        elt = 9;
        ++n;
        return true;
      }
      if(elt == 8) {
        elt = 10;
        ++n;
        return true;
      }
      return false;
    });
    EXPECT_EQ(0,nRemoved);
    EXPECT_EQ(2,n);
  }
  
  
  {
    int n = 0;
    int nRemoved = a.forEach([&n](int & elt) {
      if(elt == 10) {
        ++n;
        return false; // delete the first element
      }
      return true;
    });
    EXPECT_EQ(1,nRemoved);
    EXPECT_EQ(1,n);
  }
  
  {
    int n = 0;
    int nRemoved = a.forEach([&n](int & elt) {
      if(elt == 9) {
        ++n;
        return true;
      }
      return false;
    });
    EXPECT_EQ(0,nRemoved);
    EXPECT_EQ(1,n);
  }
  
  EXPECT_TRUE(a.tryInsert(3)); // add another element
  
  {
    int n = 0;
    int nRemoved = a.forEach([&n](int & elt) {
      if(elt == 9) {
        ++n;
        return true;
      }
      if(elt == 3) {
        ++n;
        return true;
      }
      return false;
    });
    EXPECT_EQ(0,nRemoved);
    EXPECT_EQ(2,n);
  }
  {
    int n = 0;
    int nRemoved = a.forEach([&n](int & elt) {
      ++n;
      return false;
    });
    EXPECT_EQ(2,nRemoved);
    EXPECT_EQ(2,n);
  }
}

namespace imajuscule {
  
  // To make the test more relevant,
  // we use this type whose reads / writes are not atomic
  struct RepeatInt {
    RepeatInt() {
      for(size_t i=0; i<is.size(); ++i) {
        is[i] = i;
      }
    }

    RepeatInt(int i) {
      set(i);
    }
    
    int get() const {
      return is[0];
    }
    
    void set(int i) {
      is.fill(i);
    }
    
    bool isValid() const {
      for(auto v : is) {
        if(v!= is[0]) {
          return false;
        }
      }
      return true;
    }
    
    bool operator < (RepeatInt const & other) const {
      return is[0] < other.is[0];
    }
    
  private:
    std::array<int,1000> is;
  };
}

template<int arraySz, typename T>
void testLFASlowConsumer() {
  using namespace imajuscule;
  
  LG(INFO,"using array size %d", arraySz);
  
  lockfree::scmp::static_vector<T> a(arraySz);
  
  constexpr auto nthreads = 20;
  
  std::vector<std::thread> v_threads;
  v_threads.reserve(nthreads);
  
  constexpr auto nNumbers = 10000;
  std::vector<int> values(nNumbers);
  std::iota(values.begin(), values.end(), 0);
  
  for(int i=0; i<nthreads; ++i) {
    // These are the producers : they add values to the array in a random order.
    v_threads.emplace_back(std::thread{[&a,&values] () {
      auto v = values;
      Shuffle(v);
      LG(INFO,"X");
      
      for(auto e : v) {
        bool res;
        do {
          res = a.tryInsert({e});
        } while(!res); // repeat until we succeed.
      }
      LG(INFO,".");
    }});
  }
  
  // I am the consumer
  std::vector<int> received;
  std::vector<bool> valid;
  
  received.reserve(nNumbers * nthreads);
  valid.reserve(nNumbers * nthreads);
  
  int nDeleted = 0;
  
  while(received.size() < received.capacity()) {
    nDeleted += a.forEach([&received, &valid](auto & ri) {
      
      if(rand() % 100) {
        return true; // 99% of the time, we skip the value
      }
      if constexpr (std::is_same_v<T, RepeatInt>) {
        valid.push_back(ri.isValid());
        received.push_back(ri.get());
      }
      else {
        valid.push_back(true);
        received.push_back(ri);
      }
      return false; // delete the element
    });
  }
  
  EXPECT_EQ(nDeleted, received.size());
  
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
void testLFAFastConsumer() {
  using namespace imajuscule;
  
  LG(INFO,"using array size %d and scheduling", arraySz);
  
  lockfree::scmp::static_vector<T> a(arraySz);
  
  constexpr auto nthreads = 20;
  
  std::vector<std::thread> v_threads;
  v_threads.reserve(nthreads);
  
  constexpr auto nNumbers = 10000;
  std::vector<int> values(nNumbers);
  std::iota(values.begin(), values.end(), 0);
  std::atomic<bool> run = false;
  std::atomic<bool> run2 = true;
  
  for(int i=0; i<nthreads; ++i) {
    // These are the producers : they add values to the array in a random order.
    v_threads.emplace_back(std::thread{[&a,&values,&run] () {
      auto v = values;
      Shuffle(v);
      LG(INFO,"X");
      
      for(auto e : v) {
        bool res;
        do {
          if(run) {
            res = a.tryInsert({e});
          }
          else {
            res = false;
          }
        } while(!res); // repeat until we succeed.
      }
      LG(INFO,".");
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
  
  int nDeleted = 0;

  while(received.size() < received.capacity()) {
    nDeleted += a.forEach([&received, &valid](auto & ri) {
      if constexpr (std::is_same_v<T, RepeatInt>) {
        valid.push_back(ri.isValid());
        received.push_back(ri.get());
      }
      else {
        valid.push_back(true);
        received.push_back(ri);
      }
      return false; // delete the element
    });
  }
  
  EXPECT_EQ(nDeleted, received.size());
  
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
void testLFA() {
  using namespace imajuscule;
  testLFASlowConsumer<arraySz, int>();
  testLFAFastConsumer<arraySz, int>();
  
  testLFASlowConsumer<arraySz, RepeatInt>();
  testLFAFastConsumer<arraySz, RepeatInt>();
}

TEST(StaticVectorLockFreeSCMP, multiThread) {
  testLFA<1>();
  testLFA<10>();
  testLFA<100>();
  testLFA<1000>();
  testLFA<10000>();
  testLFA<100000>();
}

namespace imajuscule {
  struct Destructible {

    Destructible() {
      countLiveObjects()++;
    }
    
    Destructible ( const Destructible & ) {
      countLiveObjects()++;
    }

    ~Destructible() {
      countLiveObjects()--;
    }
    
    static int & countLiveObjects() {
      static int n = 0;
      return n;
    }
  };
  
  static inline int nLiveObjects() {
    return Destructible::countLiveObjects();
  }
}
TEST(StaticVectorLockFreeSCMP, OnRemovalAssignFromDefault_UP) {
  using namespace imajuscule;
  
  EXPECT_EQ(0,nLiveObjects());
  {
    lockfree::scmp::static_vector<std::unique_ptr<Destructible>, lockfree::scmp::OnRemoval::AssignFromDefault> a(10);
    
    a.tryInsert(std::make_unique<Destructible>());
    EXPECT_EQ(1,nLiveObjects());
    a.tryInsert(std::make_unique<Destructible>());
    EXPECT_EQ(2,nLiveObjects());
    
    a.forEach([](auto &){ return false; });
    EXPECT_EQ(0,nLiveObjects());
  }
  EXPECT_EQ(0,nLiveObjects());
}

TEST(StaticVectorLockFreeSCMP, OnRemovalAssignFromDefault) {
  using namespace imajuscule;

  EXPECT_EQ(0,nLiveObjects());
  {
    lockfree::scmp::static_vector<Destructible, lockfree::scmp::OnRemoval::AssignFromDefault> a(10);
    
    a.tryInsert({});
    EXPECT_EQ(10,nLiveObjects());
    a.tryInsert({});
    EXPECT_EQ(10,nLiveObjects());
    
    // remove all elements
    a.forEach([](auto){ return false; });
    EXPECT_EQ(10,nLiveObjects());
  }
  EXPECT_EQ(0,nLiveObjects());
}


TEST(StaticVectorLockFreeSCMP, OnRemovalDoNothing_UP) {
  using namespace imajuscule;
  
  EXPECT_EQ(0,nLiveObjects());
  {
    lockfree::scmp::static_vector<std::unique_ptr<Destructible>, lockfree::scmp::OnRemoval::DoNothing> a(10);
    
    a.tryInsert(std::make_unique<Destructible>());
    EXPECT_EQ(1,nLiveObjects());
    a.tryInsert(std::make_unique<Destructible>());
    EXPECT_EQ(2,nLiveObjects());
    
    a.forEach([](auto &){ return false; });
    EXPECT_EQ(2,nLiveObjects()); // the removal hasn't changed the number of live objects
    
    // the next 2 inserts will succeed and not change the number of live objects
    EXPECT_TRUE(a.tryInsert(std::make_unique<Destructible>()));
    EXPECT_EQ(2,nLiveObjects());
    EXPECT_TRUE(a.tryInsert(std::make_unique<Destructible>()));
    EXPECT_EQ(2,nLiveObjects());

    // subsequent inserts will change the number of live objects
    EXPECT_TRUE(a.tryInsert(std::make_unique<Destructible>()));
    EXPECT_EQ(3,nLiveObjects());
    EXPECT_TRUE(a.tryInsert(std::make_unique<Destructible>()));
    EXPECT_EQ(4,nLiveObjects());
  }
  EXPECT_EQ(0,nLiveObjects());
}

TEST(StaticVectorLockFreeSCMP, OnRemovalDoNothing) {
  using namespace imajuscule;
  
  EXPECT_EQ(0,nLiveObjects());
  {
    lockfree::scmp::static_vector<Destructible, lockfree::scmp::OnRemoval::DoNothing> a(10);
    
    a.tryInsert({});
    EXPECT_EQ(10,nLiveObjects());
    a.tryInsert({});
    EXPECT_EQ(10,nLiveObjects());
    
    // remove all elements
    a.forEach([](auto){ return false; });
    EXPECT_EQ(10,nLiveObjects());

    a.tryInsert({});
    EXPECT_EQ(10,nLiveObjects());
    a.tryInsert({});
    EXPECT_EQ(10,nLiveObjects());
  }
  EXPECT_EQ(0,nLiveObjects());
}


