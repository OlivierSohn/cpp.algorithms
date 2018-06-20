
namespace imajuscule {
  
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
  
}
