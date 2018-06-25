

namespace imajuscule::lockfree::scmp {

  /*

   'fifo' is a lock-free single consumer, multiple producer fifo queue
   of fixed-capacity, optimized for the uncontended case (for example,
   we don't take any precautions regarding false sharing).

 -- Thread-safety --

   Multiple threads are allowed to call 'tryEnqueue'.

   A single thread at a time is allowed to call 'tryDequeue'.

 -- Lifecycle of the elements of the underlying container --

   When a 'fifo' is constructed,
    the elements of the underlying container are default-constructed.

   When a new element is queued,
    the corresponding element of the underlying container is either
     'assigned to' or 'move-assigned to' from the new element.

   When an element is dequeued,
    depending on the 'OnRemoval' template parameter,
     the corresponding element of the underlying container is either:
      left as-is, or
      'assigned to' from a default-constructed element.

   When a 'fifo' is destructed,
    the elements of the underlying container are destructed.

   */

/*
     'fifo' synopsis

namespace imajuscule::lockfree::scmp {

template<typename T, OnRemoval DP=OnRemoval::DoNothing>
struct fifo {

   // types:
   using value_type = T;

   // construction:
   fifo(int capacity);

   // thread-safe element enqueue, fails of the queue is full:
   bool tryEnqueue(value_type v);

   // (not thread-safe) element dequeue, returns an empty object if the queue is empty:
   Optional<value_type> tryDequeue();
};

} // imajuscule::lockfree::scmp

*/



  template<typename T, OnRemoval DP=OnRemoval::DoNothing>
  struct fifo {
    using value_type = T;

    static_assert(DP == OnRemoval::DoNothing ||
                  (std::is_constructible_v<value_type> && std::is_assignable_v<value_type&,value_type>));

    static constexpr auto maxCapacity = std::numeric_limits<uint16_t>::max() - 1 - 1;

    // @ param capacity : The number of elements in the queue when it is full.
    //                    This is not counting the "separator" element.
    //                    If you pass a capacity that is bigger than 'maxCapacity'
    //                      it will be silently trimmed to maxCapacity.
    //                    If you pass a capacity that is smaller than 1
    //                      it will be silently trimmed to 1.
    fifo(int capacity) :
    s(trimCapacity(capacity+1)), // we add one because when the queue is full, it has 's.size()-1' full slots.
    nextReadnextWrite(0)
    {
      for(auto it = states(), end = states_end(); it != end; ++it) {
        it->store(SlotState::Empty, std::memory_order_relaxed);
      }
    }

    // This method can be called concurrently from many threads.
    bool tryEnqueue(value_type v) {
      uint32_t rw = nextReadnextWrite.load(std::memory_order_acquire);
      while(1) {
        nonatomicReadWrite cur = getReadWrite(rw);
        if(isFull(cur)) {
          return false;
        }
        nonatomicReadWrite next = advanceWrite(cur);
        if(!nextReadnextWrite.compare_exchange_strong(rw,ungetReadWrite(next), std::memory_order_acq_rel)) {
          // failure : retry (note that rw now holds the new value)
          continue;
        }
        // success
        auto writeIndex = cur.second;
        if constexpr (std::is_move_assignable_v<value_type>) {
          values()[writeIndex] = std::move(v);
        }
        else {
          values()[writeIndex] = v;
        }
        states()[writeIndex].store(SlotState::Full, std::memory_order_release);
        break;
      }
      return true;
    }

    /*
     Thread safety:
     This method should be called from a single thread at a time.

     Applies the function f to every dequeued element.

     Returns the number of elements that were dequeued.
     */
    template<typename F>
    int dequeueAll(F f) {
      int res = 0;
      while(tryDequeue(f)) {
        ++res;
      }
      return res;
    }

    /*
     Thread safety:
     This method should be called from a single thread at a time.

     Applies the function f to the dequeued element.

     Returns the number of elements that were dequeued.
     */
    template<typename F>
    bool tryDequeue(F f) {
      uint32_t rw = nextReadnextWrite.load(std::memory_order_acquire);
      nonatomicReadWrite cur = getReadWrite(rw);
      if(isEmpty(cur)) {
        return false;
      }

      // Since this method can only be called from a single thread at a time,
      // we are able to dequeue an element from position 'cur.first' (unless the writer
      // didn't finish writing the element).

      {
        auto readIndex = cur.first;
        auto & st = states()[readIndex];
        auto currentState = SlotState::Full;
        if(!st.compare_exchange_strong(currentState, SlotState::Empty, std::memory_order_acq_rel)) {
          // the writer didn't finish writing the element.
          return false;
        }
        auto & val = values()[readIndex];
        f(val);

        if constexpr (DP == OnRemoval::AssignFromDefault) {
          val = {};
        }
      }

      nonatomicReadWrite next = advanceRead(cur);
      uint32_t nextRw = ungetReadWrite(next);
      while(!nextReadnextWrite.compare_exchange_strong(rw,nextRw, std::memory_order_acq_rel)) {
        // failure : retry

        // rw now holds the new value for reads and writes.
        // The read value cannot have changed
        // because this method is called by a single thread at a time.

        // set the 'write' bits to the new value.
        nextRw = (rw & 0xFFFF0000) | static_cast<uint32_t>(next.first);
      }
      // success
      return true;
    }


  private:
    enum class SlotState : int {
      Full      =   0b0, // the slot is full.
      Empty     = 0b100, // the slot is empty.
    };

    static_assert(std::atomic<SlotState>::is_always_lock_free);

    // We use 'DistantPairArray<...>' instead of 'std::vector<std::pair<...>>'
    // because it has a better memory layout.
    using slots = DistantPairArray<std::atomic<SlotState>,T>;

    auto * states() { return s.firsts(); }
    auto * states_end() { return s.firsts_end(); }
    auto * values() { return s.seconds(); }
    auto * values_end() { return s.seconds_end(); }

    using nonatomicReadWrite = std::pair<uint16_t,uint16_t>;

    static constexpr size_t trimCapacity(int capacity) {
      if(capacity < 1) {
        capacity = 1; // forbid queues of 0 element
      }
      if(capacity > std::numeric_limits<uint16_t>::max() - 1) {
        capacity = std::numeric_limits<uint16_t>::max() - 1; // -1 because we sometimes add 1 and compare with s.size()
      }
      return capacity;
    }

    // We could keep indexes on different cachelines to avoid false sharing.
    // But since this fifo is used for cases where the queue is almost always empty,
    // it is not obvious that it would be better on average.

    // this atomic holds the read index in its low bits , and the write index in its high bits.
    std::atomic<uint32_t> nextReadnextWrite;

    slots s;

    // empty case:
    //
    //     nextWrite
    //     v
    // --------------------------
    //     ^
    //     nextRead
    //
    // full case:
    //
    //    nextWrite
    //    v
    // --------------------------
    //     ^
    //     nextRead


    bool isFull(nonatomicReadWrite rw) {
      auto afterWrite = rw.second+1;
      return (rw.first == 0 && afterWrite==s.size()) || (rw.first == afterWrite);
    }

    static bool isEmpty(nonatomicReadWrite rw) {
      return rw.first == rw.second;
    }

    nonatomicReadWrite advanceWrite(nonatomicReadWrite rw) {
      auto afterWrite = rw.second+1;
      if(afterWrite == s.size()) {
        afterWrite = 0;
      }
      return {rw.first,afterWrite};
    }

    nonatomicReadWrite advanceRead(nonatomicReadWrite rw) {
      auto afterRead = rw.first+1;
      if(afterRead == s.size()) {
        afterRead = 0;
      }
      return {afterRead,rw.second};
    }

    static nonatomicReadWrite getReadWrite(uint32_t rw) {
      return {
        static_cast<uint16_t>(rw),
        static_cast<uint16_t>(rw >> 16)
      };
    }

    static uint32_t ungetReadWrite(nonatomicReadWrite rw) {
      return rw.first + (static_cast<uint32_t>(rw.second) << 16);
    }

  };

}
