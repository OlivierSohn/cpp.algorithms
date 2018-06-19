namespace imajuscule {

  struct LockCtrl {
    LockCtrl( std::atomic_bool & l ) noexcept : l(l) {}

    void lock() noexcept {
      bool current = false;
      while (!l.compare_exchange_weak(current,true,std::memory_order_acq_rel)) {
        std::this_thread::yield();
      }
    }

    void unlock() noexcept {
      l.store(false, std::memory_order_release);
    }
  private:
    std::atomic_bool & l;
  };


  class LockGuard {
  public:
    LockGuard( std::atomic_bool & l ) noexcept : ctrl(l) {
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
