namespace imajuscule {

  struct LockCtrl {
  public:
    LockCtrl(std::atomic_flag & f) noexcept : lock_(f) {}
    
    LockCtrl(const LockCtrl &) = delete;
    LockCtrl &operator=(const LockCtrl &) = delete;
    
    void lock() noexcept {
      while (lock_.test_and_set(std::memory_order_acquire)) {
        std::this_thread::yield();
      }
    }
    
    void unlock() noexcept { lock_.clear(std::memory_order_release); }
    
  private:
    std::atomic_flag & lock_;
  };


  class LockGuard {
  public:
    LockGuard( std::atomic_flag & l ) noexcept : ctrl(l) {
      ctrl.lock();
    }
    ~LockGuard() noexcept {
      ctrl.unlock();
    }
  private:
    LockCtrl ctrl;

    LockGuard(const LockGuard &) = delete;
    LockGuard & operator = (const LockGuard &) = delete;
  };
}
