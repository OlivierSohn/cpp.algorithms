

namespace imajuscule {
  namespace lockfree {
    
    enum class SlotState : int {
      Full      =   0b0, // the slot is full.
      SoonFull  =   0b1, // the slot is being filled.
      SoonEmpty =  0b10, // the slot is being emptied.
      Empty     = 0b100, // the slot is empty.
    };

    // We use 'PairArray<...>' instead of 'std::vector<std::pair<...>>'
    // for its more optimal memory layout.
    template<typename T>
    using slots = PairArray<std::atomic<SlotState>,T>;

    template<typename T>
    auto * states(slots<T> & s) { return s.firsts(); }
    template<typename T>
    auto * states_end(slots<T> & s) { return s.firsts_end(); }
    template<typename T>
    auto * values(slots<T> & s) { return s.seconds(); }
    template<typename T>
    auto * values_end(slots<T> & s) { return s.seconds_end(); }

    /* If the single consumer thread is more performance-critical
     than the producer threads, and if your code doesn't rely on
     destructors being called promptly, use DestructionPolicy::Lazy. */
    enum class DestructionPolicy {
      // The object destructor is called when its slot is emptied
      // (i.e in the single consumer thread).
      Prompt,
      
      // The object destructor is called when its slot is reused
      // (i.e in a producer thread), or when the whole array is destructed.
      Deferred,
    };
    
      
    /*
     
     "Single consumer - multiple producer" lock-free array, where
     elements are not required to have atomic read / write operations.
     
     The memory overhead compared to std::array is sizeof(std::atomic<int>)
     per array element.
     
     !!! As explained in its documentation, the 'forEach' method should be called
     from a single consumer thread, else data races will occur.

     */
    template<typename T, DestructionPolicy DP=DestructionPolicy::Deferred>
    struct Array {

      using value_type = T;

      static_assert(DP == DestructionPolicy::Deferred || std::is_trivially_constructible_v<value_type>);

      Array(int size) :
        s(size)
      {
        for(auto it = states(s), end = states_end(s); it != end; ++it) {
          it->store(SlotState::Empty, std::memory_order_release);
        }
        upperBoundCount.store(0,std::memory_order_release);
      }
      
      // Called by producer threads.
      //
      // Returns true if the value was added.
      bool addValue(value_type v) {
        auto begin = states(s);
        auto end = states_end(s);
        
        bool retry = false;
        do {
          for(auto state = begin; state!=end; ++state) {
            auto current = SlotState::Empty;
            if(state->compare_exchange_strong(current, SlotState::SoonFull)) {
              // update the corresponding value
              auto index = state - begin;
              values(s)[index] = v;
              
              // and change the atomic flag (by design, no other thread can have changed it in-between).
              Assert(*state == SlotState::SoonFull);
              
              // Using memory_order_seq_cst so that the counter increment will be seen before
              // the Full state, and hence upperBoundCount will actually be an upper bound.
              upperBoundCount.fetch_add(1, std::memory_order_seq_cst);
              state->store(SlotState::Full, std::memory_order_seq_cst);
              return true;
            }
            // there is no 'SoonEmpty' state when DP == DestructionPolicy::Deferred
            if constexpr (DP == DestructionPolicy::Prompt) {
              if(current == SlotState::SoonEmpty) {
                retry = true;
              }
            }
          }
        } while(retry);
        return false;
      }
      
      // !!! Producer threads are not allowed to call this function.
      //
      // This method should be called only by the single consumer.
      //
      // Elements present in the array ** at the beginning of the call ** are
      //   guaranteed to be traversed.
      // Elements added concurrently using 'addValue' may or may not be traversed.
      //
      // When the function passed as argument returns false,
      //   the value is removed from the array.
      template<typename F>
      void forEach(F f) {
        int nMax = upperBoundCount.load(std::memory_order_seq_cst);
        if(nMax == 0) {
          // the array is empty
          return;
        }
        Assert(nMax > 0);
        // the array may be non empty.
        auto begin = states(s);
        auto end = states_end(s);
        auto begin_values = values(s);
        auto val = begin_values;
        int nSeen = 0;
        int nRemoved = 0;
        for(auto state = begin; state!=end; ++state, ++val) {
          if(state->load(std::memory_order_seq_cst) != SlotState::Full) {
            continue;
          }
          ++nSeen;

          if(!f(*val)) {
            ++nRemoved;
            // there is no 'SoonEmpty' state when DP == DestructionPolicy::Deferred
            if constexpr (DP == DestructionPolicy::Prompt) {
              state->store(SlotState::SoonEmpty, std::memory_order_release);
              *val = {};
            }
            state->store(SlotState::Empty, std::memory_order_release);
          }
          
          if(nSeen == nMax) {
            // check again, as producers may have added elements.
            nMax = upperBoundCount.load(std::memory_order_seq_cst);
            if(nSeen == nMax) {
              break;
            }
            // Some elements have been added by producers
            // between the beginning of the iteration and now,
            // so we continue. Note that we may need to traverse
            // the whole array, in case the added element(s) have
            // been missed by this loop.
          }
        }

        // relaxed here because only the consumer cares about ordering.
        upperBoundCount.fetch_add(-nRemoved, std::memory_order_relaxed);
      }
      

    private:
      slots<value_type> s;

      // An upper bound of the count of full slots.
      std::atomic<int> upperBoundCount;
    };
  }
}
