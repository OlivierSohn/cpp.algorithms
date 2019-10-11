namespace imajuscule
{
  /*
   Not thread safe.
   
   Allows to control the variation of a value, to avoid brusk changes.
   */
  template<typename T>
  struct smoothVar {
    smoothVar(T v) {
      smoothly(v,0);
    }
    
    /*
     * Sets the target value and the count of steps of the interpolation.
     */
    void smoothly(T v, int nSteps) {
      if(nSteps == 0) {
        current = target = v;
        increment = 0;
        return;
      }
      if(unlikely(nSteps < 0)) {
        nSteps = -nSteps;
      }
      target = v;
      increment = std::abs(target - current) / static_cast<T>(nSteps);
      Assert(increment >= 0);
    }
    
    /*
     * Gets the current value, using linear interpolation to perform the transition.
     */
    T step() {
      Assert(increment >= 0);
      if(current == target) {
      }
      else if(current < target-increment) {
        current += increment;
      }
      else if(current > target+increment) {
        current -= increment;
      }
      else {
        current = target;
      }
      return current;
    }
      
      T getWithoutStepping() const {
          return current;
      }
  private:
    T current, target, increment;
  };
  
}
