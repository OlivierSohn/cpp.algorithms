
/*
 When using std::array<pair<A,B>> where A and B have different alignment constraints,
 padding can be used in the std::pair, and memory will be wasted.
 
 'PairArray<A,B>' fixes this issue by not using std::pair, and uses a single
 contiguous memory buffer containing:
   * first the As or Bs, depending on which has a bigger alignment constraint,
   * then the others (Bs or As).
 */

namespace imajuscule {

  enum class AssignSeconds;
  
  namespace detail {

  enum class Order {
    As_Bs, // As are first
    Bs_As  // Bs are first
  };
  
  enum class ArrayLocality {
      Local
    , Distant
  };
  
  template<typename A, typename B, ArrayLocality Loc, int szLocalArray>
  struct PairArray {
    static_assert(szLocalArray >= 0);
    static_assert(szLocalArray == 0 || Loc == ArrayLocality::Local);

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

    struct DistantArray {

      DistantArray(int countPairs) :
      count(countPairs),
      data(allocateBuffer(countPairs))
      {}
      
      ~DistantArray() {
        detail::deallocate_aligned_memory(data);
      }
      
      BeginT * buf() { return data; }
      BeginT const * buf() const { return data; }
      auto size() const { return count; }

    private:
      int count;
      BeginT * data;
    };
    
    struct LocalArray {
      BeginT * buf() { return reinterpret_cast<BeginT*>(a.data()); }
      BeginT const * buf() const { return reinterpret_cast<BeginT const*>(a.data()); }
      constexpr auto size() const { return szLocalArray; }
    private:

      alignas(buf_align) std::array<char,bufferSize(szLocalArray)> a;
    };
    static_assert(alignof(LocalArray) == buf_align);

    using Arr = std::conditional_t<Loc==ArrayLocality::Distant, DistantArray, LocalArray>;
    
  public:
    
    template<ArrayLocality L = Loc
    , std::enable_if_t<
        L == ArrayLocality::Distant &&
        std::is_constructible_v<A> &&
        std::is_constructible_v<B>>...>
    PairArray(int countPairs) : arr{countPairs}
    {
      static_assert(Loc == ArrayLocality::Distant);
      initialize();
    }
    
    template<typename A_, typename B_, ArrayLocality L = Loc
    , std::enable_if_t<L == ArrayLocality::Distant>...>
    PairArray(int countPairs, A_ && a, B_ && b) : PairArray(countPairs)
    {
      fill(a,b);
    }
    
    template<ArrayLocality L = Loc
    , std::enable_if_t<
        L == ArrayLocality::Local &&
        std::is_constructible_v<A> &&
        std::is_constructible_v<B>>...>
    PairArray()
    {
      static_assert(Loc == ArrayLocality::Local);
      initialize();
    }
    
    template<typename A_, typename B_, ArrayLocality L = Loc
    , std::enable_if_t<L == ArrayLocality::Local>...>
    PairArray( A_ && a, B_ && b) : PairArray()
    {
      fill(a,b);
    }
    
    template<typename A_, typename B_, ArrayLocality L = Loc
    , std::enable_if_t<L == ArrayLocality::Local>...>
    PairArray(A_ && a, std::array<B_, szLocalArray> & bs)
    {
      initialize_a(std::move(a));
      initialize_bs(bs);
    }

    template<typename A_, typename B_, ArrayLocality L = Loc
    , std::enable_if_t<L == ArrayLocality::Local>...>
    PairArray(std::array<A_, szLocalArray> & as, B_ && b)
    {
      initialize_as(as);
      initialize_b(std::move(b));
    }
    
    template<typename A_, typename B_, ArrayLocality L = Loc
    , std::enable_if_t<L == ArrayLocality::Local>...>
    PairArray(std::array<A_, szLocalArray> & as, std::array<B_, szLocalArray> & bs)
    {
      initialize_as(as);
      initialize_bs(bs);
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
    
  private:
    void initialize()
    {
      for(auto it = lefts(), end = lefts_end(); it != end; ++it) {
        new (it) BeginT;
      }
      for(auto it = rights(), end = rights_end(); it != end; ++it) {
        new (it) EndT;
      }
    }

    template<typename A_>
    void initialize_as(std::array<A_, szLocalArray> & as)
    {
      auto a = as.begin();
      for(auto it = firsts(), end = firsts_end(); it != end; ++it, ++a) {
        new (it) A{*a};
      }
    }
    
    template<typename B_>
    void initialize_bs(std::array<B_, szLocalArray> & bs)
    {
      auto b = bs.begin();
      for(auto it = seconds(), end = seconds_end(); it != end; ++it, ++b) {
        new (it) B{*b};
      }
    }
    
    template<typename A_>
    void initialize_a(A_ && a)
    {
      for(auto it = firsts(), end = firsts_end(); it != end; ++it) {
        new (it) A;
        *it = a;
      }
    }
    
    template<typename B_>
    void initialize_b(B_ && b)
    {
      for(auto it = seconds(), end = seconds_end(); it != end; ++it) {
        new (it) B;
        *it = b;
      }
    }
    
  public:
    ~PairArray() {
      for(auto it = lefts(), end = lefts_end(); it != end; ++it) {
        it->~BeginT();
      }
      for(auto it = rights(), end = rights_end(); it != end; ++it) {
        it->~EndT();
      }
    }


    int size() const {
      return arr.size();
    }
    
    // should we provide an indexing operator that works with std::pair<A,B> ?
    //    T & operator[](size_t n) { return ...; };
    //    auto operator[](size_t n) const -> std::pair<A,B> { return ...; };
    
    // returns the pointer to the portion of the buffer containing the first elements
    A * firsts() {
      if constexpr (order == Order::As_Bs) {
        return arr.buf();
      }
      else {
        return reinterpret_cast<A*>(arr.buf() + arr.size());
      }
    }
    
    A const * firsts() const {
      if constexpr (order == Order::As_Bs) {
        return arr.buf();
      }
      else {
        return reinterpret_cast<A*>(arr.buf() + arr.size());
      }
    }
    
    // returns the pointer to the portion of the buffer containing the second elements
    B * seconds() {
      if constexpr (order == Order::Bs_As) {
        return arr.buf();
      }
      else {
        return reinterpret_cast<B*>(arr.buf() + arr.size());
      }
    }
    
    B const * seconds() const {
      if constexpr (order == Order::Bs_As) {
        return arr.buf();
      }
      else {
        return reinterpret_cast<B*>(arr.buf() + arr.size());
      }
    }
    
    // returns the pointer to the end of the portion of the buffer containing the first elements
    A * firsts_end() {
      if constexpr (order == Order::As_Bs) {
        return arr.buf() + arr.size();
      }
      else {
        return reinterpret_cast<A*>(reinterpret_cast<char*>(arr.buf()) + arr.size() * sizeof_AB);
      }
    }
    
    // returns the pointer to the end of the portion of the buffer containing the second elements
    B * seconds_end() {
      if constexpr (order == Order::Bs_As) {
        return arr.buf() + arr.size();
      }
      else {
        return reinterpret_cast<B*>(reinterpret_cast<char*>(arr.buf()) + arr.size() * sizeof_AB);
      }
    }
    
    A& corresponding(B const &b) {
      auto i = &b - seconds();
      return firsts()[i];
    }
    
    B& corresponding(A const &a) {
      auto i = &a - firsts();
      return seconds()[i];
    }

  private:
    Arr arr;
    
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
  } // detail
  
  // The array memory is local
  template<typename A, typename B, int N>
  using LocalPairArray = detail::PairArray<A,B,detail::ArrayLocality::Local,N>;

  // The array memory is distant (a pointer is used to point to the array)
  template<typename A, typename B>
  using DistantPairArray = detail::PairArray<A,B,detail::ArrayLocality::Distant, 0>;
  
  // range iteration
  template <typename T>
  class RangeIter
  {
  public:
    RangeIter(T* collection, size_t size) :
    mCollection(collection), mSize(size)
    {
    }
    
    T* begin() { return &mCollection[0]; }
    T* end() { return &mCollection[mSize]; }
    
  private:
    T* mCollection;
    size_t mSize;
  };

  template<typename T>
  auto firsts(T && a) {
    return RangeIter(a.firsts(), a.size());
  }
  
  template<typename T>
  auto seconds(T && a) {
    return RangeIter(a.seconds(), a.size());
  }

}
