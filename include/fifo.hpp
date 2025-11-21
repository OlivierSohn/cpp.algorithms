/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule {
    /*
     * Non-thread-safe, FIFO queue of unlimited size, based on a vector, which performs no
     * dynamic memory allocation / deallocation within its methods: the
     * trick is that the caller of methods whose names are 'unsafe'-prefixed
     * are responsible for giving the queue its new underlying container,
     * if the queue is full.
     *
     * This separation of concern allows to use the queue in a real-time context, where
     * dynamic allocations are forbidden.
     *
     * @see safeTryEmplace
     *
     * The object T should be default-constructible.
    */
    template< typename T >
    struct fifo {
      using container = std::vector<T>;
      using value_type = typename container::value_type;
      using iterator = typename container::iterator;

      static constexpr auto min_capacity = 1;

      fifo(int initialCapacity = 10) {
        // we must have at least one element in the vector so that full() has a valid implementation.
        int capacity = 1 + std::max(0, initialCapacity);
        v.resize(capacity);
        read = write = v.begin();
      }

      /*
       * It is very important that this function doesn't free or allocate any memory.
       */
      void reset() {
        if(read <= write) {
          for (auto i=read; i!=write; ++i) {
            *i = {};
          }
        }
        else {
          for (auto i=read, end = v.end(); i!=end; ++i) {
            *i = {};
          }
          for (auto i=v.begin(); i!=write; ++i) {
            *i = {};
          }
        }
        read = write = v.begin();
      }

      T & front() {
          return *read;
      }

      T const & front() const {
          return *read;
      }

      T & back() {
          if(write == v.begin()) {
              return v.back();
          }
          return *(write-1);
      }

      T const & back() const {
          if(write == v.begin()) {
              return v.back();
          }
          return *(write-1);
      }

      bool empty() const {
          return read == write;
      }

      bool full() const {
          return size() == v.size()-1;
      }

      /*
       * Returns the additional capacity that should be added to the current underlying container
       * in order to support adding 'nAdds' element to the queue.
       *
       * 0 is returned when the queue can hold the elements without growing, else
       * a strictly positive number is returned.
       */
      int shouldGrow(int nAdds) const {
        int nElems = size();
        int maxCount = v.size() - 1;
        int diff = nElems + nAdds - maxCount;
        if( diff > 0 ) {
          return diff;
        }
        return 0;
      }

      std::size_t size() const {
          if(read <= write) {
              return std::distance(read, write);
          }
          return v.size() - std::distance(write, read);
      }

      int32_t underlyingContainerCapacity() const { return v.capacity(); }

      /*
       * Before calling this method, please make sure the fifo is not full.
       *
       * @See safeTryEmplace.
       */
      void emplace( T&& value ) {
        Assert(!full());
        *write = std::move(value);
        inc(write);
      }

      /*
       * Removes the last added element.
       */
      void cancel_emplace() {
        Assert(!empty());
        dec(write);
        *write = {}; // to destroy the object
      }

      void pop() {
        Assert(!empty());
        *read = {}; // to destroy the object
        inc(read);
      }

      void swap( fifo& other ) noexcept {
        std::swap(v, other.v);
        std::swap(read, other.read);
        std::swap(write, other.write);
      }

      /*
       * This method will swap the underlying container if the passed container has
       * a strictly bigger capacity than the original one.
       *
       * @Returns true if the swap has been done, false otherwise.
       *
       * If the swap was done, the original storage is returned in c.
       */
      bool trySwapUnderlyingContainer(container & c) {
        Assert(c.size() == 0);

        if(c.capacity() <= v.capacity()) {
          return false;
        }

        // move the elements to the new container
        if(read <= write) {
          std::move(read,write,std::back_insert_iterator<container>(c));
        }
        else {
          std::move(read,v.end(),std::back_insert_iterator<container>(c));
          std::move(v.begin(),write,std::back_insert_iterator<container>(c));
        }
        write = c.end();
        c.resize(c.capacity());
        Assert(write < c.end()); // because c.capacity() > v.capacity()
        read = c.begin();
        std::swap(c,v);
        return true;
      }

    private:
      iterator read, write;

      container v;

      void inc (iterator&i) {
        ++i;
        if(i==v.end()) {
          i = v.begin();
        }
      }

      void dec (iterator&i) {
        if(i==v.begin()) {
          i = v.end();
        }
        --i;
      }
    };

  /*
   * Adds an element to the fifo queue (using emplace) after having made sure
   * that the fifo queue has enough room for it.
   *
   * The 'CanRealloc' template parameter controls whether the underlying container may be reallocated
   * or not (if we are lockfree, we don't reallocate)
   *
   * @param l : a potential lock (no-op if we are lockfree) protecting the queue writes.
   *        The allocation / deallocation of memory will happen outside the lock scope,
   *        which makes this function suitable for use in a real-time context, provided that
   *        the lock can configure the thread priorities to avoid priority inversion.
   *
   * @returns whether or not the element was emplaced
   */
  template<CanRealloc canRealloc, typename T, typename Lock>
  [[nodiscard]] bool safeTryEmplace(fifo<T> & q, Lock & l, T && v) {
    bool res = reserveAndLock<canRealloc>(1,q,l);
    if(res) {
      q.emplace(std::move(v));
    }

    l.unlock();

    return res;
  }

} // NS imajuscule
