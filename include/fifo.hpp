/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule {
    /*
     * FIFO queue of unlimited size, based on a vector, which performs no
     * dynamic memory allocation / deallocation within its methods: the
     * trick is that the caller of methods whose names are 'unsafe'-prefixed
     * are responsible for giving the queue its new underlying container,
     * if the queue is full.
     *
     * This separation of concern allows to use the queue in a real-time context, where
     * dynamic allocations are forbidden.
     *
     * @see safeEmplace
     *
     * The object T should be default-constructible.
    */
    template< typename T >
    struct fifo {
      using container = std::vector<T>;
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
      
      std::size_t underlying_container_capacity() const { return v.capacity(); }
      
      /*
       * Before calling this method, please make sure the fifo is not full.
       *
       * @See safeEmplace.
       */
      void emplace( T&& value ) {
        Assert(!full());
        *write = std::move(value);
        inc(write);
      }

      void pop() {
        Assert(!empty());
        *read = {}; // to destroy the object
        inc(read);
      }

      void swap( fifo& other )
        noexcept (
          std::is_nothrow_swappable<container>::value &&
          std::is_nothrow_swappable<iterator>::value
                 ) {
        std::swap(v, other.v);
        std::swap(read, other.read);
        std::swap(write, other.write);
      }
      
      /*
       * The original storage will be returned in c.
       *
       * If the vector passed as parameter has not a bigger capacity than the current vector, the function does nothing.
       */
      void swapStorage(container & c) {

        Assert(c.capacity() > v.capacity());
        Assert(c.size() == 0);
        
        if(c.capacity() <= v.capacity()) {
          return;
        }
        
        // move the elements to the new container
        if(read <= write) {
          std::move(read,write,std::back_insert_iterator(c));
        }
        else {
          std::move(read,v.end(),std::back_insert_iterator(c));
          std::move(v.begin(),write,std::back_insert_iterator(c));
        }
        write = c.end();
        c.resize(c.capacity());
        Assert(write < c.end()); // because c.capacity() > v.capacity()
        read = c.begin();
        std::swap(c,v);
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
    };
  
  namespace detail {
    /*
     * In a multiple-producers / multiple consumers context, the non-const
     * queue methods will be protected by a lock.
     */
    template<typename T, typename Lock>
    void makeSomeRoom(int n, std::size_t cap1, fifo<T> & q, Lock & l)
    {
      {
        static constexpr auto growth_factor = 2;
        std::vector<T> newStorage;
        
        // allocation happens outside the lock scope:
        newStorage.reserve(std::min(cap1 + 1, growth_factor * cap1));
        
        l.lock();
        
        auto cap2 = q.underlying_container_capacity();
        
        if(newStorage.capacity() > cap2) {
          q.swapStorage(newStorage);
        }

        l.unlock();

        // deallocation happens outside the lock scope:
        newStorage = std::vector<T>{}; // we make vector desctruction explicit for clarity but removing this line is just as good.
      }
      
      l.lock();

      if(auto amount = q.shouldGrow(n)) {
        auto cap3 = q.underlying_container_capacity();
        l.unlock();
        
        makeSomeRoom(n, cap3+amount, q, l);
      }
    }
  }

  /*
   * Reallocates, if needed, the underlying container of the queue so that 'n'
   * additional elements can be pushed. When this function returns, the lock 'l'
   * is being taken, hence it will be safe to add the elements right after this call,
   * and then the caller should release the lock.
   *
   * 'n' : the number of elements that we want to be able to add to the queue.
   */
  template<typename T, typename Lock>
  void reserveAndLock(int n, fifo<T> & q, Lock & l) {
    l.lock();
    
    if(auto amount = q.shouldGrow(n)) {
      auto cap1 = q.underlying_container_capacity();
      
      l.unlock();
      
      detail::makeSomeRoom(n, cap1+amount, q, l);
      // the lock was taken by the previous call.
    }
  }

  /*
   * Adds an element to the fifo queue (using emplace) after having made sure
   * that the fifo queue has enough room for it.
   *
   * The potential lock protecting the queue writes is represented
   * by the two functions passed as argument.
   *
   * If needed, the allocation / deallocation of memory will happen outside the lock scope,
   * which makes this function suitable for use in a real-time context, provided that
   * lock / unlock functions can configure the thread priorities to avoid priority inversion.
   */
  template<typename T, typename Lock>
  void safeEmplace(fifo<T> & q, Lock & l, T && v) {
    reserveAndLock(1,q,l);
    
    q.emplace(std::move(v));

    l.unlock();
  }

} // NS imajuscule
