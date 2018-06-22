
namespace imajuscule {
  
  ////////////////////////////////////////////////
  
  struct BigAlignment {
    BigAlignment() = default;
    BigAlignment(double d) {
      v = d;
    }
    alignas (64) double v;
  };
  
  struct SmallAlignment {
    double v;
  };
  
  static_assert(alignof(SmallAlignment)<alignof(BigAlignment));
  static_assert(alignof(SmallAlignment)>=alignof(void*));
  
  
  //////////////////////////////////////////////////
  // A type whose reads / writes are not atomic
  struct RepeatInt {
    RepeatInt() {
      for(size_t i=0; i<is.size(); ++i) {
        is[i] = i;
      }
    }
    
    RepeatInt(int i) {
      set(i);
    }
    
    int get() const {
      return is[0];
    }
    
    void set(int i) {
      is.fill(i);
    }
    
    bool isValid() const {
      for(auto v : is) {
        if(v!= is[0]) {
          return false;
        }
      }
      return true;
    }
    
    bool operator < (RepeatInt const & other) const {
      return is[0] < other.is[0];
    }
    
  private:
    std::array<int,1000> is;
  };
  
  /////////////////////////////////////////
  struct Destructible {
    
    Destructible() {
      countLiveObjects()++;
    }
    
    Destructible ( const Destructible & ) {
      countLiveObjects()++;
    }
    
    ~Destructible() {
      countLiveObjects()--;
    }
    
    static int & countLiveObjects() {
      static int n = 0;
      return n;
    }
  };
  
  static inline int nLiveObjects() {
    return Destructible::countLiveObjects();
  }
  
  ////////////////////////////////////////////
}
