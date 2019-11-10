namespace imajuscule
{
template<int polyDegree, typename T>
T poly(T v) {
    static_assert(polyDegree >= 0);
    static_assert(polyDegree <= 6);
    if constexpr (polyDegree == 0) {
        return 1;
    }
    else if constexpr(polyDegree == 1) {
        return v;
    }
    else if constexpr(polyDegree == 2) {
        return v*v;
    }
    else if constexpr(polyDegree == 3) {
        return v*v*v;
    }
    else if constexpr(polyDegree == 4) {
        auto v2 = v*v;
        return v2*v2;
    }
    else if constexpr(polyDegree == 5) {
        auto v2 = v*v;
        return v2*v2*v;
    }
    else if constexpr(polyDegree == 6) {
        auto v2 = v*v;
        return v2*v2*v2;
    }
}
  /*
   Not thread safe.
   
   Allows to control the variation of a value, to avoid brusk changes.
   */
  template<typename T, int polyDegree>
  struct smoothVar {
      static_assert(polyDegree >= 1);
      
    smoothVar(T v)
      : current(v)
      , target(v)
      , increment(0)
      {}
    
    /*
     * Sets the target value and the count of steps of the interpolation.
     */
    void smoothly(T v, int nSteps) {
      if(nSteps == 0) {
        increment = 0;
        target = v;
        current.setX(v);
        return;
      }
      target = v;
      increment = std::abs(target - current.getX()) / static_cast<T>(std::abs(nSteps));
      Assert(increment >= 0);
    }
    
    /*
     * Gets the current value, using linear interpolation to perform the transition.
     */
    T step() {
      Assert(increment >= 0);
        auto curX = current.getX();
      if(curX == target) {
      }
      else if(curX < target-increment) {
        curX += increment;
      }
      else if(curX > target+increment) {
        curX -= increment;
      }
      else {
        curX = target;
      }
      current.setX(curX);
      return getWithoutStepping();
    }
      
      T getWithoutStepping() const {
          return current.getY();
      }
  private:
      struct PrivateValue {
          explicit PrivateValue(T x) :
          x(x) {
              computeY();
          }
          void setX(T val) {
              if(x == val) {
                  return;
              }
              x = val;
              computeY();
          }
          
          T getX() const {
              return x;
          }
          
          T getY() const {
              return y;
          }
      private:
          T x;
          T y;
          
          void computeY() {
              y = poly<polyDegree>(x);
          }
      };
      PrivateValue current;
    T target, increment;
  };
  
}
