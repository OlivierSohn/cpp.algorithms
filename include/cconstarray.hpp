
namespace imajuscule
{
  template<typename T>
  struct CConstArray {
    using value_type = T;

    std::size_t size() const { return sz; }
    CConstArray(T const * buf, int sz_) : buf(buf), sz(sz_) {
      if(sz < 0) {
        Assert(0);
        sz = 0;
      }
    }

    auto begin() const { return buf; }
    auto end() const { return buf + sz; }

    auto const & operator [] (int i) const {
      Assert(i >= 0);
      Assert(i < sz);
      return buf[i];
    }

  private:
    T const *buf;
    int sz;
  };

  template<int Sz, typename T>
  struct CArray {
    using value_type = T;
    
    constexpr std::size_t size() const { return Sz; }
    CArray(T * buf) : buf(buf) {
    }
    
    auto begin()     { return buf; }
    auto end()       { return buf + Sz; }
    auto begin() const { return buf; }
    auto end() const { return buf + Sz; }
    
    auto const & operator [] (int i) const {
      Assert(i >= 0);
      Assert(i < Sz);
      return buf[i];
    }
    
  private:
    T *buf;
  };

}
