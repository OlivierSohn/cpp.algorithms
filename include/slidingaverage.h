
namespace imajuscule
{
  // TODO use advanced traits for efficient parameter passing
  
// also privodes sliding max / min values
  template <class T>
  struct slidingAverage
  {
    using ParameterType = T;
    
    slidingAverage(slidingAverage&&) = default;
    slidingAverage& operator=(slidingAverage&&) = default;
    
    slidingAverage(size_t size) :
    values(size)
    {}
    
    void feed(ParameterType val) {
      if(firstFeed) {
        firstFeed = false;
        values.feedAndFill(val);
      }
      else {
        values.feed(val);
      }
    }
    
    T compute() const {
      return std::accumulate(values.begin(), values.end(), T{}) / safe_cast<float>(values.size());
    }
    T computeMax() const {
      return *std::max_element(values.begin(), values.end());
    }
    
    void reset() {
      values.reset();
      firstFeed = true;
    }
    
    auto size() const {
      return values.size();
    }
    
  private:
    cyclic<T> values;
    bool firstFeed = true;
  };
  
  
  /*
   * the window function parameter is between 0 (out of the window) and 1 (middle of the window)
   * the 0 value is never used
   */
  template <class T>
  struct slidingWindowedAverage
  {
    using ParameterType = T;
    
    slidingWindowedAverage(slidingWindowedAverage&&) = default;
    slidingWindowedAverage& operator=(slidingWindowedAverage&&) = default;
    
    template<typename WINDOW>
    slidingWindowedAverage(size_t size, WINDOW f) :
    values(size),
    window(size)
    {
      assert(std::abs(f(1.f) - 1.f) < 0.0001f);
      
      float n_steps;
      float mid_idx;
      if(size%2) {
        // even
        n_steps = 1.f + (size-1)/2;
        mid_idx = n_steps;
      }
      else {
        // odd
        n_steps = .5f + size/2;
        mid_idx = n_steps;
      }
      
      auto it = window.begin();
      auto end = window.end();
      auto step = 1;
      for(; it!=end; ++it, ++step) {
        float real_step = step;
        if(real_step > mid_idx) {
          real_step = 2.f * mid_idx - real_step;
        }
        
        *it = f(real_step / n_steps);
        assert(*it <= 1.f);
      }
    }
    
    void feed(ParameterType val) {
      if(firstFeed) {
        firstFeed = false;
        values.feedAndFill(val);
      }
      else {
        values.feed(val);
      }
    }
    
    T compute() const {
      T sum{};
      auto it = window.begin();
      values.for_each([&](auto v) {
        sum += v * *it;
        ++it;
      });
      return sum / safe_cast<float>(values.size());
    }
    
    void reset() {
      values.reset();
      firstFeed = true;
    }
    
    auto size() const {
      return values.size();
    }
    
    auto getWindowAverage() const {
      return std::accumulate(window.begin(), window.end(), T{}) / window.size();
    }
  private:
    cyclic<T> values;
    std::vector<T> window;
    bool firstFeed = true;
  };
}
