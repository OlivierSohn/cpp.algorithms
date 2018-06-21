

namespace imajuscule {

  enum class Atomicity {
    Yes,
    No
  };
  
  template<Atomicity A, typename T>
  struct maybeAtomic;
  
  template<typename T>
  struct maybeAtomic<Atomicity::Yes, T> {
    using UnderlyingType = T;
    using type = std::atomic<T>;
    static constexpr auto atomic = Atomicity::Yes;

    static_assert(type::is_always_lock_free);

    static bool compareExchangeStrong(std::atomic<T> & value, T from, T to, std::memory_order order) {
      T fr = from;
      return value.compare_exchange_strong(fr,to,order);
    }
    
    [[nodiscard]] static T read(std::atomic<T> const & value, std::memory_order order) {
      return value.load(order);
    }
    
    static void write(std::atomic<T> & value, T v, std::memory_order order) {
      value.store(v,order);
    }
  };
  
  template<typename T>
  struct maybeAtomic<Atomicity::No, T> {
    using UnderlyingType = T;
    using type = T;
    static constexpr auto atomic = Atomicity::No;

    static bool compareExchangeStrong(T & value, T from, T to, std::memory_order) {
      if(value == from) {
        value = to;
        return true;
      }
      else {
        return false;
      }
    }

    [[nodiscard]] static T read(T value, std::memory_order) {
      return value;
    }
    
    static void write(T & value, T v, std::memory_order) {
      value = v;
    }
  };
  
  
  template<typename U, typename Enable = void>
  struct FlagTraits;
  
  template<typename T>
  struct FlagTraits<std::atomic<T>, std::enable_if_t<std::atomic<T>::is_always_lock_free> > {
    using traits = maybeAtomic<Atomicity::Yes, std::atomic<T>>;
  };
  
  template<typename T>
  struct FlagTraits<T, std::enable_if_t<std::is_trivially_copyable_v<T>>> {
    using traits = maybeAtomic<Atomicity::No, T>;
  };
  
  


  
}
