

namespace imajuscule::lockfree::scmp {

    /*
     
     'static_vector' is a lock-free variable-size container with fixed-capacity.
     
   -- Thread-safety --
     
     Multiple threads are allowed to call 'tryInsert' to add new elements.
     
     A single thread is allowed to call 'forEach' to read / modify / remove
     the elements.
     
   -- Lifecycle of the elements of the underlying container --
     
     When a 'static_vector' is constructed,
      the elements of the underlying container are default-constructed.

     When a new element is added to a 'static_vector', using 'tryInsert',
      the corresponding element of the underlying container is either
       'assigned to' or 'move-assigned to' from the new element.

     When an element is removed from a 'static_vector', using 'forEach',
      depending on the 'OnRemoval' template parameter,
       the corresponding element of the underlying container is either:
        left as-is, or
        'assigned to' from a default-constructed element.

     When a 'static_vector' is destructed,
      the elements of the underlying container are destructed.

   -- Memory usage --
     
     The memory overhead compared to std::vector of the same capacity
     is 'sizeof(std::atomic<int>)' per element.
     
     */
    
/*
     'static_vector' synopsis

 namespace imajuscule::lockfree::scmp {

template<typename T, OnRemoval DP=OnRemoval::DoNothing>
struct static_vector {
 
     // types:
     using value_type = T;

     // construction:
     static_vector(int capacity);
     
     // thread-safe element insertion:
     bool tryInsert(value_type v);

     // (not thread-safe) element read/write access and removal:
     template<typename F>
     int forEach(F f);
};

} // imajuscule::lockfree::scmp

*/

    enum class OnRemoval {
      // The object is assigned from a default-constructed element.
      AssignFromDefault,
      
      // The object is left as-is.
      DoNothing,
    };

    template<typename T, OnRemoval DP=OnRemoval::DoNothing>
    struct static_vector {
      using value_type = T;

      static_assert(DP == OnRemoval::DoNothing ||
                    (std::is_constructible_v<value_type> && std::is_assignable_v<value_type&,value_type>));

      static_vector(int capacity) :
        s(capacity)
      {
        for(auto it = states(), end = states_end(); it != end; ++it) {
          it->store(SlotState::Empty, std::memory_order_release);
        }
        upperBoundCount.store(0,std::memory_order_release);
      }
      
      // Many threads can call this method concurrently.
      //
      // Returns true if the insertion succeeded, false otherwise.
      bool tryInsert(value_type v) {
        auto begin = states();
        auto end = states_end();

        bool retry = false;
        do {
          for(auto state = begin; state!=end; ++state) {
            auto current = SlotState::Empty;
            if(state->compare_exchange_strong(current, SlotState::SoonFull)) {
              // update the corresponding value
              auto index = state - begin;
              if constexpr (std::is_move_assignable_v<value_type>) {
                values()[index] = std::move(v);
              }
              else {
                values()[index] = v;
              }
              
              // and change the atomic flag (by design, no other thread can have changed it in-between).
              Assert(*state == SlotState::SoonFull);
              
              // Using memory_order_seq_cst so that the counter increment will be seen before
              // the Full state, and hence upperBoundCount will actually be an upper bound.
              upperBoundCount.fetch_add(1, std::memory_order_seq_cst);
              state->store(SlotState::Full, std::memory_order_seq_cst);
              return true;
            }
            // there is no 'SoonEmpty' state when DP == OnRemoval::DoNothing
            if constexpr (DP == OnRemoval::AssignFromDefault) {
              if(current == SlotState::SoonEmpty) {
                // we don't wait, but we will retry if the insertion didn't happen in a subsequent slot
                retry = true;
              }
            }
          }
        } while(retry);
        return false;
      }

      // This method allows to read, modify, remove the elements present in the static_vector
      //   "at the beginning of the call".
      //
      // A single thread (at a time) is allowed to call this method. Reentrancy is not supported.
      //
      // Elements that are potentially added to the static_vector "during the call"
      // are not guaranteed to be traversed, but will be traversed during subsequent calls.
      //
      // @param f: A callable that will be called once per element. The element
      //           can be passed:
      //             - by value (a copy will be made)
      //             - or by const reference (no copy will be made)
      //             - or by reference (it is then possible to modify the element in the callable)
      //           If the callable returns true, the element is kept in the static_vector, else
      //           it is removed.
      //
      // @returns the number of elements that were removed from the static_vector.
      template<typename F>
      int forEach(F f) {
        int nMax = upperBoundCount.load(std::memory_order_seq_cst);
        if(nMax == 0) {
          // the static_vector is empty
          return 0;
        }
        Assert(nMax > 0);
        // the static_vector may be non empty.
        auto begin = states();
        auto end = states_end();
        auto begin_values = values();
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
            // there is no 'SoonEmpty' state when DP == OnRemoval::DoNothing
            if constexpr (DP == OnRemoval::AssignFromDefault) {
              state->store(SlotState::SoonEmpty, std::memory_order_release);
              *val = {};
            }
            state->store(SlotState::Empty, std::memory_order_release);
          }
          
          if(nSeen == nMax) {
            nMax = upperBoundCount.load(std::memory_order_seq_cst);
            if(nSeen == nMax) {
              break;
            }
            // Some elements have been added
            // between the beginning of the iteration and now,
            // so we continue. Note that we may need to traverse
            // the whole underlying container, in case the added element() have
            // been missed by this loop.
          }
        }

        // relaxed here because only the consumer cares about ordering.
        upperBoundCount.fetch_add(-nRemoved, std::memory_order_relaxed);
        return nRemoved;
      }
      
    private:
      enum class SlotState : int {
        Full      =   0b0, // the slot is full.
        SoonFull  =   0b1, // the slot is being filled.
        SoonEmpty =  0b10, // the slot is being emptied.
        Empty     = 0b100, // the slot is empty.
      };
      
      static_assert(std::atomic<SlotState>::is_always_lock_free);
      
      // We use 'DistantPairArray<...>' instead of 'std::vector<std::pair<...>>'
      // because it has a better memory layout.
      using slots = DistantPairArray<std::atomic<SlotState>,T>;
      
      // An upper bound of the count of full slots.
      std::atomic<int> upperBoundCount;
      slots s;
    
      auto * states() { return s.firsts(); }
      auto * states_end() { return s.firsts_end(); }
      auto * values() { return s.seconds(); }
      auto * values_end() { return s.seconds_end(); }
    };
      
}
