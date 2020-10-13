#define IMJ_DEBUG_COMPRESSOR 0

namespace imajuscule::audio {

  template<typename T>
  struct SafeSince {
      SafeSince() {
          onUnsafe();
      }

      void onUnsafe() {
          since = 0;
          max_safe_amplitude = 0.;
      }

      void onSafe(T amplitude) {
          ++since;
          max_safe_amplitude = std::max(amplitude, max_safe_amplitude);
      }

      int getSafeSince() const { return since; }
      T getSafeMaxAmplitude() const { return max_safe_amplitude; }
  private:
      int since;
      T max_safe_amplitude;
  };

  /* This compressor adapts the signal volume to its amplitude.
   *
   * To avoid having having compression level changes eventhough the signal
   * didn't change much, we use hysteresis.
   *
   * if the multiplication factor is different from the target,
   * it changes at most by a factor of compressMult or uncompressMult per frame.
   * But if the compressed signal would clip, the multiplication factor is instantaneously changed
   * to the smallest value that avoids clipping.
   */
  template<typename T>
  struct Compressor {
    // do not expand the signal, only compress
    static constexpr T maxTargetMultiplicator = 1.;
    
    static constexpr T limit = 0.999;
    static constexpr T high = 0.7;
    static constexpr T low  = 0.5;

    static constexpr T compressMult   = 0.997;

    // = 2^(1/22050) so that the volume doubles every half second
    static constexpr T uncompressMult = 1.00003; // TODO take sample rate into account

    // The signal needs to be low for at least that many frames
    // to be uncompressed
    static constexpr int safeDuration = 100000; // TODO take sample rate into account

    Compressor() {
#if IMJ_DEBUG_COMPRESSOR
        std::cout << "Compressor targetMultiplicator = " << targetMultiplicator
                  << " (initial value)" << std::endl;
        std::cout << "Compressor multiplicator = " << multiplicator
                  << " (initial value)" << std::endl;
#endif
    }

    template<typename S>
    void feed(S & signal) {
      static_assert(std::is_same_v<typename S::value_type, T>);

      Assert(targetMultiplicator <= 1.);
      Assert(targetMultiplicator > 0.);

      T rawAmplitude{};
      for(auto s : signal) {
        rawAmplitude = std::max(rawAmplitude, std::abs(s));
      }

      // Adjust the _target_ compression level,
      // taking only into account the amplitude of the signal "when the compression level will have reached its target",
      if(rawAmplitude * targetMultiplicator > high) {
        // the signal is too loud.
        // we adjust the target such that the rawAmplitude would correspond to medium level:
        //   targetMultiplicator * rawAmplitude = medium
        //   targetMultiplicator                = medium / rawAmplitude
        Assert(rawAmplitude);
        targetMultiplicator = std::min(maxTargetMultiplicator, medium / rawAmplitude);
#if IMJ_DEBUG_COMPRESSOR
        std::cout << "Compressor targetMultiplicator = " << targetMultiplicator
                  << " (signal was too high)" << std::endl;
#endif
        safeSince.onUnsafe();
      }
      else if(rawAmplitude * targetMultiplicator < low) {
        // the signal is low
        safeSince.onSafe(rawAmplitude);
        if(safeSince.getSafeSince() > safeDuration) {
          // The signal has been too low for a long time
          if (targetMultiplicator < maxTargetMultiplicator) {
            T const maxRawAmplitude = safeSince.getSafeMaxAmplitude();
            // we adjust the target such that the maxRawAmplitude would correspond to medium level:
            //   targetMultiplicator * maxRawAmplitude = medium
            //   targetMultiplicator                   = medium / maxRawAmplitude
            if (maxRawAmplitude) {
              targetMultiplicator = std::min(maxTargetMultiplicator, medium / maxRawAmplitude);
            } else {
              targetMultiplicator = maxTargetMultiplicator;
            }
#if IMJ_DEBUG_COMPRESSOR
            std::cout << "Compressor targetMultiplicator = " << targetMultiplicator
                      << " (signal low for " << safeDuration << " frames : targetAmplitude = " << rawAmplitude * targetMultiplicator << " < " << low << ")" << std::endl;
#endif
          }
          safeSince.onUnsafe();
        }
      }
      else {
        // the signal is medium, keep the current compression level
        safeSince.onUnsafe();
      }

      // slowly adjust the compression level to match targetMultiplicator
      if(likely(multiplicator == targetMultiplicator)) {
        // do nothing
      }
      else {
        if(multiplicator < targetMultiplicator) {
          // TODO(OS) we could use another strategy to have fixed duration release (to avoid cases where the duration is too long)
          // (for attack, we use the clipping strategy + fixed multiplicator so there is no need to have fixed duration
          //  because the duration will never be too long)
          multiplicator = std::min(multiplicator * uncompressMult, targetMultiplicator);
#if IMJ_DEBUG_COMPRESSOR
          std::cout << "Adjusting compressor multiplicator = " << multiplicator
                    << " (lower than target)" << std::endl;
#endif
        }
        else {
          multiplicator = std::max(multiplicator * compressMult, targetMultiplicator);
#if IMJ_DEBUG_COMPRESSOR
          std::cout << "Adjusting compressor multiplicator = " << multiplicator
                    << " (bigger than target)" << std::endl;
#endif
        }
      }

      // immediately adjust the compression level if the compressed signal clips
      T amplitude = rawAmplitude * multiplicator;
      if(amplitude > limit) {
        multiplicator = limit / rawAmplitude;
#if IMJ_DEBUG_COMPRESSOR
        std::cout << "Forcing compressor multiplicator = " << multiplicator
                  << " (to avoid clipping : amplitude = " << amplitude << " > " << limit << ")" << std::endl;
#endif
        safeSince.onUnsafe();
      }

      // compress the signal.
      for(auto & s : signal) {
        s *= multiplicator;
      }
    }

    // needed for white box tests
    T getCompressionLevel() const { return multiplicator; }
    T getTargetCompressionLevel() const { return targetMultiplicator; }
    int getSafeSince() const { return safeSince.getSafeSince(); }
  private:
    T multiplicator = 1.;
    T targetMultiplicator = 1.;
    SafeSince<T> safeSince;
    const T medium = std::exp( (std::log(high) + std::log(low)) / 2.);
  };
}
