


template<typename SL, typename USED>
void testLock(USED & used) {
    
    constexpr auto nloops = 1000000;
    constexpr auto nthreads = 2;
    
    unsigned long long shared = 0;
    
    std::vector<std::thread> v_threads;
    v_threads.reserve(nthreads);
    
    for(int i=0; i<nthreads; ++i) {
        v_threads.emplace_back(std::thread{[&shared, &used] () {
            for(int i=0; i<nloops; ++i) {
                SL l(used);
                ++shared;
            }
        }});
    }
    
    for(auto & t: v_threads) {
        t.join();
    }
    
    EXPECT_EQ(nthreads * nloops, shared);
}

// todo compare performances in different contexts
TEST(Locks, compare) {
  using namespace imajuscule;
  {
    LG(INFO,"mutex");
    std::mutex m;
    testLock<std::lock_guard<std::mutex>>(m);
  }
  {
    LG(INFO,"atomic_flag");
    std::atomic_flag f = ATOMIC_FLAG_INIT;
    testLock<LockGuard>(f);
  }
  
}
