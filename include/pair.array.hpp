
/*
 When using std::array<pair<A,B>> where A and B have different alignment constraints,
 padding can be used in the std::pair, and memory will be wasted.
 
 'PairArray<A,B>' fixes this issue by not using std::pair, and uses a single
 contiguous memory buffer containing:
   * first the As or Bs, depending on which has a bigger alignment constraint,
   * then the others (Bs or As).
 */

namespace imajuscule {
  
  
  enum class Order {
    As_Bs, // As are first
    Bs_As  // Bs are first
  };

  
  template<typename A, typename B>
  struct PairArray {

    // We store the elements that have the biggest alignment constraint first.
    static constexpr auto order =
      alignof(A) >= alignof(B) ?
        Order::As_Bs :
        Order::Bs_As;

  private:
    static_assert(is_power_of_two(alignof(A)));
    static_assert(is_power_of_two(alignof(B)));

    static constexpr auto buf_align = std::max(sizeof(void*),std::max(alignof(A), alignof(B)));
    static constexpr auto sizeof_AB = sizeof(A) + sizeof(B);
    
    // the type of elements that are stored at the beginning of the array
    using BeginT = std::conditional_t<order == Order::As_Bs, A, B>;
    // the type of elements that are stored at the end of the array
    using EndT = std::conditional_t<order == Order::As_Bs, B, A>;
    
    // Assumes that there will be no gap between elements in the buffer:
    // this is true because the smaller alignment is divisible by the bigger one,
    // and elements having the bigger alignment are stored first in the buffer.
    static constexpr size_t bufferSize(int countPairs) {
      return countPairs * sizeof_AB;
    }
    static BeginT * allocateBuffer(int countPairs) {
      return reinterpret_cast<BeginT*>(detail::allocate_aligned_memory(buf_align, bufferSize(countPairs)));
    }
    
  public:
    // TODO should I call constructors? destructors? (on every element)

    PairArray(int countPairs) :
    count(countPairs)
    , buf(allocateBuffer(countPairs))
    {
      for(auto it = lefts(), end = lefts_end(); it != end; ++it) {
        new (it) BeginT;
      }
      for(auto it = rights(), end = rights_end(); it != end; ++it) {
        new (it) EndT;
      }
    }
    
    ~PairArray() {
      for(auto it = lefts(), end = lefts_end(); it != end; ++it) {
        it->~BeginT();
      }
      for(auto it = rights(), end = rights_end(); it != end; ++it) {
        it->~EndT();
      }
      detail::deallocate_aligned_memory(buf);
    }

    template<typename A_, typename B_>
    PairArray(int countPairs, A_ && a, B_ && b) : PairArray(countPairs)
    {
      fill(a,b);
    }
    
    template<typename A_, typename B_>
    void fill(A_ && a, B_ && b)
    {
      for(auto it = firsts(), end = firsts_end(); it != end; ++it) {
        *it = a;
      }
      for(auto it = seconds(), end = seconds_end(); it != end; ++it) {
        *it = b;
      }
    }
        
    int size() const {
      return count;
    }
    
    // should we provide an indexing operator that works with std::pair<A,B> ?
    //    T & operator[](size_t n) { return ...; };
    //    auto operator[](size_t n) const -> std::pair<A,B> { return ...; };
    
    // returns the pointer to the portion of the buffer containing the first elements
    A * firsts() {
      if constexpr (order == Order::As_Bs) {
        return buf;
      }
      else {
        return reinterpret_cast<A*>(buf + count);
      }
    }
    
    // returns the pointer to the portion of the buffer containing the second elements
    B * seconds() {
      if constexpr (order == Order::Bs_As) {
        return buf;
      }
      else {
        return reinterpret_cast<B*>(buf + count);
      }
    }
    
    // returns the pointer to the end of the portion of the buffer containing the first elements
    A * firsts_end() {
      if constexpr (order == Order::As_Bs) {
        return buf + count;
      }
      else {
        return reinterpret_cast<A*>(reinterpret_cast<char*>(buf) + count * sizeof_AB);
      }
    }
    
    // returns the pointer to the end of the portion of the buffer containing the second elements
    B * seconds_end() {
      if constexpr (order == Order::Bs_As) {
        return buf + count;
      }
      else {
        return reinterpret_cast<B*>(reinterpret_cast<char*>(buf) + count * sizeof_AB);
      }
    }
    
  private:
    // allocate_aligned_memory
    int count;
    BeginT * buf;
    
    auto * lefts() {
      if constexpr (order == Order::As_Bs) {
        return firsts();
      }
      else {
        return seconds();
      }
    }

    auto * rights() {
      if constexpr (order == Order::As_Bs) {
        return seconds();
      }
      else {
        return firsts();
      }
    }

    auto * lefts_end() {
      if constexpr (order == Order::As_Bs) {
        return firsts_end();
      }
      else {
        return seconds_end();
      }
    }
    
    auto * rights_end() {
      if constexpr (order == Order::As_Bs) {
        return seconds_end();
      }
      else {
        return firsts_end();
      }
    }
  };
}
