

namespace imajuscule {

  enum class Synchronization {
    // a single thread can access the elements of the container, modify them, and remove them.
    // multiple threads can add elements to the container
    Lockfree_SingleConsumerMultipleProducer,

    // every operation on the container should happen from a single thread at a time.
    SingleThread
  };
  
  template<Synchronization S, typename ...Args>
  struct static_vector_s;
  
  template<typename ...Args>
  struct static_vector_s<Synchronization::SingleThread, Args...> {
    using type = singlethread::static_vector<Args...>;
  };

  template<typename ...Args>
  struct static_vector_s<Synchronization::Lockfree_SingleConsumerMultipleProducer, Args...> {
    using type = lockfree::scmp::static_vector<Args...>;
  };

  template<Synchronization S, typename ...Args>
  using static_vector = typename static_vector_s<S, Args...>::type;

}
