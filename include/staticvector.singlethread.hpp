

namespace imajuscule {
  namespace singlethread {
    
    /*
     
     'static_vector' is a variable-size container with fixed-capacity.
     
     A synchronized, lock-free version exists in 'staticvector.lockfree.hpp'.
     
   -- Thread-safety --
     
     A single thread is allowed to call 'tryInsert' to add new elements.
     
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

namespace imajuscule {
 
template<typename T, OnRemoval DP=OnRemoval::DoNothing>
struct static_vector {
 
     // types:
     using value_type = T;

     // construction:
     static_vector(int capacity);
     
     // (not thread-safe) element insertion:
     bool tryInsert(value_type v);

     // (not thread-safe) element read/write access and removal:
     template<typename F>
     int forEach(F f);
};

} // imajuscule

*/

    enum class OnRemoval {
      // The object is assigned from a default-constructed element.
      AssignFromDefault,
      
      // The object is left as-is.
      DoNothing,
    };

    namespace detail {
      enum class SlotState : uint8_t {
        Full      =   0b0, // the slot is full.
        Empty     = 0b100, // the slot is empty.
      };

      // We use 'PairArray<...>' instead of 'std::vector<std::pair<...>>'
      // because it has a better memory layout.
      template<typename T>
      using slots = PairArray<SlotState,T>;

      template<typename T>
      auto * states(slots<T> & s) { return s.firsts(); }
      template<typename T>
      auto * states_end(slots<T> & s) { return s.firsts_end(); }
      template<typename T>
      auto * values(slots<T> & s) { return s.seconds(); }
      template<typename T>
      auto * values_end(slots<T> & s) { return s.seconds_end(); }
    }
      
    template<typename T, OnRemoval DP=OnRemoval::DoNothing>
    struct static_vector {

      using value_type = T;

      static_assert(DP == OnRemoval::DoNothing ||
                    (std::is_constructible_v<value_type> && std::is_assignable_v<value_type&,value_type>));

      static_vector(int capacity) :
        s(capacity)
      {
        using namespace detail;
        for(auto it = states(s), end = states_end(s); it != end; ++it) {
          *it = SlotState::Empty;
        }
        count = 0;
      }
      
      // Element insertion.
      //
      // A single thread at a time is allowed to call this method.
      // Concurrent calls to 'forEach' are not supported.
      //
      // Returns true if the insertion succeeded, false otherwise.
      bool tryInsert(value_type v) {
        using namespace detail;
        auto begin = states(s);
        for(auto state = begin, end = states_end(s); state!=end; ++state) {
          if(*state == SlotState::Empty) {
            *state = SlotState::Full;
            // update the corresponding value
            auto index = state - begin;
            if constexpr (std::is_move_assignable_v<value_type>) {
              values(s)[index] = std::move(v);
            }
            else {
              values(s)[index] = v;
            }
            
            ++count;
            return true;
          }
        }
        return false;
      }
      
      // This method allows to read, modify, remove the elements present in the static_vector.
      //
      // A single thread at a time is allowed to call this method.
      // Concurrent calls to 'tryInsert' are not supported.
      // Reentrancy is not supported.
      //
      //
      // @param f: A callable that will be called once per element. The element
      //           can be passed:
      //             - by value (a copy will be made)
      //             - or by const reference (no copy will be made)
      //             - or by reference (it is then possible to modify the element in the callable)
      //           If the callable returns true, the element is kept in the static_vector, else
      //           it is removed.
      //
      // If the callable adds elements to the static_vector, they
      // are not guaranteed to be traversed, but will be traversed during subsequent calls.
      //
      // @returns the number of elements that were removed from the static_vector.
      template<typename F>
      int forEach(F f) {
        using namespace detail;
        int n = count;
        if(n == 0) {
          return 0;
        }
        Assert(n > 0);
        auto begin = states(s);
        auto end = states_end(s);
        auto begin_values = values(s);
        auto val = begin_values;
        int nSeen = 0;
        int nRemoved = 0;
        for(auto state = begin; state!=end; ++state, ++val) {
          if(*state != SlotState::Full) {
            continue;
          }
          ++nSeen;

          if(!f(*val)) // this call can modify 'count'
          {
            ++nRemoved;
            if constexpr (DP == OnRemoval::AssignFromDefault) {
              *val = {};
            }
            *state = SlotState::Empty;
          }
          
          if(nSeen == n) {
            n = count;
            if(nSeen == n) {
              break;
            }
            // Some elements have been added
            // between the beginning of the iteration and now,
            // so we continue. Note that we may need to traverse
            // the whole underlying container, in case the added element(s) have
            // been missed by this loop.
          }
        }
        
        count -= nRemoved;
        return nRemoved;
      }
      

    private:
      detail::slots<value_type> s;

      // The count of full slots.
      int count;
    };
  }
}
