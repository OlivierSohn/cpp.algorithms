

namespace imajuscule::lockfree::scmp {

  /*

   'forward_list' is a forward-linked-list with :
     - Concurrent atomic insertion at the beginning ............... (O(1))
     = Concurrent atomic removal of the last inserted element ..... (O(1))

   -- Thread-safety --

   Multiple threads call 'emplace_front', 'try_pop_front' concurrently to add / remove elements

   A single thread calls 'forEach' to traverse the list and read / modify its elements.

   -- Lifecycle of elements --

   This ensures that an element is never deleted while being traversed:

   1. Elements are constructed on dynamically allocated memory during 'emplace_front'.
   2. The next 'try_pop_front' marks the front element to be deleted.
   3. The next 'forEach' puts elements marked to be deleted in a garbage forward list.
   4. The next 'emplace_front' or 'try_pop_front' or 'collect_garbage' flushes the garbage list,
        the element is destroyed, and the dynamic memory returned to the system.

   -- See also --

   'imajuscule::lockfree::scmp::static_vector' which :
   - has an underlying container with contiguous memory,
   - is fixed-capacity,
   - may have faster element construction because the memory is statically allocated.

   */

/*
     'forward_list' synopsis

namespace imajuscule::lockfree {

// element removal:
struct Remover {
 void flagForRemoval(std::memory_order order = std::memory_order_release);
};

template<typename T>
struct forward_list {

   // types:
   using value_type = T;

   // atomic element insertion at the beginning (can be called concurrently)
   template <class... Args>
   std::pair<value_type&, Remover> emplace_front(Args&&... args);

   // atomic first element removal (can be called concurrently)
   void try_pop_front();

   // element read/write access (should be called by a single thread at a time)
   template<typename F>
   void forEach(F f);
};

 } // imajuscule::lockfree::scmp

*/
  namespace detail {
    enum class ElementFlag {
      Present,
      ShouldBeRemoved
    };

    using flag_type = std::atomic<ElementFlag>;

    template<typename T>
    struct list_elt {
      template <class... Args>
      list_elt(list_elt * elt, Args&&... args) :
      v(std::forward<Args>(args)...),
      p(elt)
      {
        flag.store(ElementFlag::Present, std::memory_order_relaxed);
      }

      flag_type flag;
      list_elt * p;
      T v;
    };

    template<typename T>
    void destroy(list_elt<T> * e) {
      auto * p = e;
      while(p) {
        auto next = p->p;
        delete p;
        p = next;
      }
    }
  }

  struct Remover {
    Remover (detail::flag_type&flag) : flag(flag) {}

    // Flag the associated element for removal.
    void flagForRemoval(std::memory_order order = std::memory_order_release) {
      using namespace detail;
      flag.store(ElementFlag::ShouldBeRemoved, order);
    }
  private:
    detail::flag_type& flag;
  };

  template<typename T>
  struct forward_list {

    using value_type = T;


  private:
    using elt = detail::list_elt<value_type>;

    struct CollectGarbageOnExit {
      CollectGarbageOnExit(forward_list&u) : u(u) {}
      ~CollectGarbageOnExit() {
        u.collect_garbage();
      }
      forward_list&u;
    };

  public:

    forward_list() {
      first.store(nullptr, std::memory_order_relaxed);
      garbage.store(nullptr, std::memory_order_relaxed);
    }

    ~forward_list() {
      clear();
    }

    // This method can be called only when the list is not traversed.
    void clear() {
      using namespace detail;
      destroy(  first.exchange(nullptr, std::memory_order_acq_rel));
      destroy(garbage.exchange(nullptr, std::memory_order_acq_rel));
    }

    // Atomic element insertion at the beginning.
    //
    // Can be called concurrently by multiple threads.
    //
    // @returns a reference to the inserted value, and a 'Remover' that
    // allows to mark the element for removal.
    template <class... Args>
    std::pair<value_type&, Remover> emplace_front(Args&&... args) {
      using namespace detail;
      CollectGarbageOnExit c(*this);

      // 1. create the new element
      auto * f = first.load(std::memory_order_acquire);
      auto * newElt = new elt(f, std::forward<Args>(args)...);
      // 2. insert it
      while(!first.compare_exchange_strong(f,newElt, std::memory_order_acq_rel)) {
        newElt->p = f;
      }
      return {
        newElt->v,
        {newElt->flag}
      };
    }

    // 'try_pop_front' removes the front element iff either:
    //   - it is the first 'try_pop_front' call since the last 'emplace_front' call
    //   - or if 'forEach' was called, and returned, since the last 'try_pop_front' call returned.
    //
    // Can be called concurrently by multiple threads.
    //
    // @returns true if an element was removed, false otherwise.
    //   Note that it can returns false on a non-empty list (cf. documentation above)
    //
    // @ Question: Why is there no 'pop_front' available?
    // @ Answer: Implementing 'pop_front' is tricky because 'forEach'
    //     concurrently re-routes the forward pointers.
    //     In my application(s), the guarantees provided by 'try_pop_front'
    //     are sufficient, so I'll take the time to write 'pop_front' only if
    //     I actually need it.
    //
    // @ Question: How can I remove a specific element?
    // @ Answer: The corresponding 'emplace_front' returns a 'Remover' : use its
    //     'flagForRemoval' method.
    bool try_pop_front() {
      using namespace detail;
      CollectGarbageOnExit c(*this);

      auto * f = first.load(std::memory_order_acquire);
      if(!f) {
        return false;
      }
      auto expected = ElementFlag::Present;
      // Note that we could use std::memory_order_release, the only downside could be that the method returns sometimes
      // true for an element that already had the flag set to 'detail::ElementFlag::ShouldBeRemoved'.
      return f->flag.compare_exchange_strong(expected, ElementFlag::ShouldBeRemoved, std::memory_order_acq_rel);
    }

    // element read/write access (should be called by a single thread at a time)
    template<typename F>
    void forEach(F func) {
      using namespace detail;
      elt * local_garbage = nullptr;
      auto * firstPtr = first.load(std::memory_order_acquire);
      elt * prevValid (nullptr);
      elt * firstValid = nullptr;
      bool need_relink = false;

      auto * f = firstPtr;
      while(f) {
        if(unlikely(f->flag.load(std::memory_order_relaxed) == ElementFlag::ShouldBeRemoved)) {
          if(local_garbage == nullptr) {
            // take ownership of 'garbage'
            local_garbage = garbage.exchange(nullptr, std::memory_order_acq_rel);
          }
          // save the pointer to the next element
          auto * next = f->p;
          // prepend f to the local_garbage
          f->p = local_garbage;
          local_garbage = f;
          // continue
          f = next;
          need_relink = true;
        }
        else {
          func(f->v); // NOTE we could allow the function to remove elements by returning false.
          if(!firstValid) {
            firstValid = f;
          }
          if(prevValid && need_relink) {
            prevValid->p = f;
          }
          prevValid = f;
          need_relink = false;
          f = f->p;
        }
      }

      if(prevValid && need_relink) {
        prevValid->p = nullptr;
      }

      if(firstValid != firstPtr) {
        f = firstPtr;
        if(!first.compare_exchange_strong(f, firstValid, std::memory_order_acq_rel)) {
          // a concurrent 'emplace_front' call happenned, which changed 'first':
          // we walk the new elements, and update the one that needs to relink to firstValid.

          while(1) {
            Assert(f); // if we don't find the element that needs to be relinked, we have a logic error.
            auto next = f->p;
            if(next == firstPtr) {
              // we found the element that needs to be relinked
              f->p = firstValid;
              break;
            }
            f = next;
          }
        }
      }

      if(local_garbage) {
        garbage.store(local_garbage, std::memory_order_release);
      }
    }

    // This method is called during 'emplace_front' and 'try_pop_front',
    // but you can also call it directly, if you have a use case where
    // that makes sense.
    //
    // This method will call 'delete', so it is not safe to use it in
    // a real-time thread.
    //
    // This method can be called concurrently.
    void collect_garbage() {
      using namespace detail;
      destroy(garbage.exchange(nullptr, std::memory_order_acq_rel));
    }

  private:
    std::atomic<elt *> first;
    std::atomic<elt *> garbage;

  };
}
