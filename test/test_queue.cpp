
TEST(RTQueue, test) {
  using Queue = atomic_queue::AtomicQueueB2<
  /* T = */ int,
  /* A = */ std::allocator<int>,
  /* MAXIMIZE_THROUGHPUT */ true,
  /* TOTAL_ORDER = */ true,
  /* SPSC = */ true
  >;
  
  for (auto capacity : std::vector<int>{100, 4096, 96000}) {
    std::cout << "min capacity:" << capacity << std::endl;

    for (int n = 1; n < 4; ++n) {
      Queue q(capacity);
      
      for (int j=0; j<10; ++j) {
        
        for (int i=0; i<q.capacity()/n; ++i) {
          ASSERT_TRUE(q.try_push(i));
        }
        
        for (int i=0; i<q.capacity()/n; ++i) {
          int j;
          ASSERT_TRUE(q.try_pop(j));
          ASSERT_EQ(i, j);
        }
        
        for (int i=0; i<q.capacity(); ++i) {
          ASSERT_TRUE(q.try_push(i));
        }
        for (int i=0; i<q.capacity(); ++i) {
          ASSERT_FALSE(q.try_push(i));
        }
        
        for (int i=0; i<q.capacity(); ++i) {
          int j;
          ASSERT_TRUE(q.try_pop(j));
          ASSERT_EQ(i, j);
        }
      }
    }
  }
}


TEST(RTQueue, threads) {
  using Queue = atomic_queue::AtomicQueueB2<
  /* T = */ int,
  /* A = */ std::allocator<int>,
  /* MAXIMIZE_THROUGHPUT */ true,
  /* TOTAL_ORDER = */ true,
  /* SPSC = */ true
  >;
  
  for (auto capacity : std::vector<int>{100, 4096, 96000}) {
    std::cout << "min capacity:" << capacity << std::endl;
    
    Queue q(capacity);
    int constexpr repetitions = 100;
    
    std::thread other([&q](){
      int n = 0;
      for (int k = 0; k < repetitions; ++k) {
        for (int i=0; i<q.capacity(); ++i) {
          int value;
          while(!q.try_pop(value)) {}
          ASSERT_TRUE(i == value);
          ++n;
          if (0 == n % 20) {
            std::this_thread::yield();
          }
          if (0 == n % (q.capacity()/10)) {
            std::this_thread::sleep_for(std::chrono::milliseconds(2));
          }
        }
      }
    });
    
    for (int k = 0; k < repetitions; ++k) {
      for (int i=0; i<q.capacity(); ++i) {
        bool res = q.try_push(i);
        if (k==0) {
          ASSERT_TRUE(res);
        }
        while(!res) {
          res = q.try_push(i);
        }
      }
    }
    
    other.join();
  }
}

