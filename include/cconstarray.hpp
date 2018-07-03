
namespace imajuscule
{
  template<typename T>
  struct CConstArray {
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

}
