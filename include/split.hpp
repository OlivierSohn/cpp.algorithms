

namespace imajuscule {

  /*
   * 'Split' stores two unsigned values in the bits of an underlying unsigned type.
   *
   * for code clarity, we let the user specify both nHigh and nLow, and
   * static_assert that "nLow + nHigh == 4*sizeof(U)".
   */
  template<typename U, U nHigh, U nLow>
  struct Split {

    static_assert(std::is_unsigned_v<U>);

    static_assert(nLow + nHigh == 8*sizeof(U));

    static_assert(nLow > 0);
    static_assert(nHigh > 0);

    static constexpr U maxHigh = std::numeric_limits<U>::max() >> nLow;
    static constexpr U maxLow  = std::numeric_limits<U>::max() >> nHigh;

    static constexpr U highMask = ~maxLow;
    static_assert(highMask == std::numeric_limits<U>::max() - maxLow);

    Split() = default;
    Split(U raw) : raw(raw) {}

    // for safety, we mask 'low' so that it won't modify the bits where
    // 'high' is stored, if it's out of bounds.
    Split(U high, U low) : raw((high << nLow) + (low & maxLow)) {
      Assert(low <= maxLow); // else, some information is lost
      Assert(high <= maxHigh); // else, some information is lost
    }

    U getHighWithZeros() const {
      return raw & highMask;
    }

    U getHigh() const {
      return raw >> nLow;
    }
    U getLow() {
      return raw & maxLow;
    }
  private:
    U raw;
  };
}
