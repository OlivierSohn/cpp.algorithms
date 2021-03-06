/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule {

    /*
    * Non-thread-safe FIFO queue containing at most a - single - element.
    *
    * It performs no dynamic allocation.
    *
    * It can be used as a replacement for std::queue when we know that there will be
    * at most one element in the queue.
    */
    template< typename T >
    struct fifo1 {
      
      void reset() {
        v = {};
      }
      
      T & front() {
        return get_value(v);
      }

      T const & front() const {
        return get_value(v);
      }

      T & back() {
        return get_value(v);
      }

      T const & back() const {
        return get_value(v);
      }

      bool empty() const { return !v; }
      bool full() const { return !empty(); }

      std::size_t size() const {
        if(v) {
          return 1;
        }
        return 0;
      }

      void push( const T& value ) {
        Assert(!v);
        v = value;
      }

      void push( T&& value ) {
        Assert(!v);
        v = std::move(value);
      }

      template< class... Args >
      decltype(auto) emplace( Args&&... args ) {
        using namespace optional_ns;
        v = make_optional<T>(std::forward<Args>(args)...);
      }
      
      /*
       * Removes the last added element.
       */
      void cancel_emplace() {
        Assert(v);
        v = {};
      }

      void pop() {
        Assert(v);
        v = {};
      }

      void swap( fifo1& other ) noexcept {
        std::swap(v, other.v);
      }

    private:
      // using an optional because: (http://en.cppreference.com/w/cpp/utility/optional)
      //
      // If an optional<T> contains a value, the value is guaranteed to be allocated
      // as part of the optional object footprint, i.e. no dynamic memory allocation
      // ever takes place.
      Optional<T> v;
    };

} // NS imajuscule
