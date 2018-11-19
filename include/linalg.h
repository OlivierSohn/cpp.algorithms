namespace imajuscule {

  template<typename T>
  struct Matrix {
    
    auto operator [] (int i) {
      return m[i];
    }

    T get(int i, int j) const {
      return v[i][j];
    }
    
    T & edit(int i, int j) {
      return v[i][j];
    }
    
    void resize(int n, int m) {
      m.resize(n);
      for(auto & v : m) {
        v.resize(m);
      }
    }
  private:
    std::vector<std::vector<T>> m;
  };
  
  template<typename T>
  auto mkIdentity(int n) {
    Matrix res;
    res.resize(n,n);
    for(int i=0; i<n; ++i) {
      res.edit(i,i) = 1;
    }
    return res;
  }
}
