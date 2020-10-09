/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule {

template <typename iterator>
struct iter_range {
    using value_type = typename std::iterator_traits<iterator>::value_type;

    iter_range(iterator begin, iterator end) : begin_(begin), end_(end) {}

    iterator begin() const {
        return begin_;
    }

    void setBegin(iterator begin) {
        begin_ = begin;
    }

    iterator end() const {
        return end_;
    }

    void setEnd(iterator end) {
        end_ = end;
    }

    int distance() const {
        return std::distance(begin_, end_);
    }

    template <typename F>
    void foreach(F && f) const {
        std::for_each(begin_, end_, f);
    }

    bool operator < (iter_range const & other) const {
        return std::distance(end_, other.begin_) >= 0;
    }

    bool follows(iter_range const & other) const {
        return begin_ == other.end_;
    }
private:
    iterator begin_, end_;
};


  // emulation of Haskell's Data.List.take, in range-based for loops
  template<typename Container>
  struct take {
    take(int n, Container & c): n(n), c(c) {}

    auto begin() const { return c.begin(); }
    auto end()   const {
      return std::min(c.end(), c.begin() + n);
    }
    bool empty() const {
      return n == 0;
    }
    std::size_t size() const {
      return n;
    }
    auto & operator [](int i) {
      return c[i];
    }
    auto const & operator [](int i) const {
      return c[i];
    }
  private:
    int n;
    Container & c;
  };

  template<typename Container>
  struct containerRange {
    using It = typename Container::iterator;
    auto begin() const { return b; }
    auto end()   const { return e; }

    containerRange(It b, It e) : b(b), e(e) {}
    Container materialize() const {
      return {b,e};
    }
    
    bool empty() const { return begin() == end(); }
    
    void setBegin(It beg) { b = beg; }
    void setEnd(It end) { e = end; }
  private:
    It b,e;
  };

  template<typename Container>
  std::pair<containerRange<Container>, containerRange<Container>> splitAt(int n, Container & c) {
    auto b = c.begin();
    auto e = c.end();
    auto m = e;
    if(n >= 0 && b+n < e) {
      m = b+n;
    }
    return {{b,m},{m,e}};
  }

  template<typename Container>
  Container withLinearFadeIn(int fadeSz, Container c) {
    // "0" and "1" values are before and after the 'fadeSz' values modified here.
    if(fadeSz < 0) {
      throw(std::logic_error("negative fade"));
    }
    if(c.size() < fadeSz) {
      throw(std::logic_error("cannot fade tiny container"));
    }
    int nIntervals = fadeSz + 1;
    using T = typename Container::value_type;
    static_assert(std::is_floating_point_v<T>);

    T step = 1 / static_cast<T>(nIntervals);
    
    for(int i=0; i<fadeSz; ++i) {
      c[i] *= step * (i+1);
    }
    
    return std::move(c);
  }
  
  
  template<typename Container>
  Container withLinearFadeOut(int fadeSz, Container c) {
    // "1" and "0" values are before and after the 'fadeSz' values modified here.
    if(fadeSz < 0) {
      throw(std::logic_error("negative fade"));
    }
    if(c.size() < fadeSz) {
      throw(std::logic_error("cannot fade tiny container"));
    }
    int nIntervals = fadeSz + 1;
    using T = typename Container::value_type;
    static_assert(std::is_floating_point_v<T>);
    
    T step = 1 / static_cast<T>(nIntervals);
    
    for(int i=0; i<fadeSz; ++i) {
      c[c.size()-1-i] *= step * (i+1);
    }

    return std::move(c);
  }
} // NS imajuscule
