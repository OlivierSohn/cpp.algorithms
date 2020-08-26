#define IMJ_DEBUG_COMPRESSOR 0

namespace imajuscule::audio {

  struct SafeSince {
      void onUnsafe() {
          since = 0;
          max_safe_amplitude = 0.f;
      }
      void onSafe(float amplitude) {
          ++since;
          max_safe_amplitude = std::max(amplitude, max_safe_amplitude);
      }

      int getSafeSince() const { return since; }
      float getSafeMaxAmplitude() const { return max_safe_amplitude; }
  private:
      int since = 0;
      float max_safe_amplitude = 0.f;
  };

  /* This compressor quickly reduces the signal volume
   * when the signal becomes too loud. and slowly restores
   * the initial volume when the signal has been below a threshold for a
   * significant period of time.
   *
   * To avoid having having compression level changes eventhough the signal
   * didn't change much, we use hysteresis.
   *
   * if the multiplication factor is different from the target,
   * it changes at most by a factor of compressMult or uncompressMult per frame.
   * But if the compressed signal would clip, the multiplication factor is instantaneously changed
   * to the smallest value that avoids clipping.
   */
  struct Compressor {
    // compress == 0.5f was making a too big change in volume.
    static constexpr float compress   = 0.8f;
    static constexpr float uncompress = 1.f / compress;

    static constexpr float limit = 0.999f;
    static constexpr float high = 0.7f;
    static constexpr float low  = 0.6f * compress;

    static constexpr float compressMult   = 0.997f;

    // = 2^(1/22050) so that the volume doubles every half second
    static constexpr float uncompressMult = 1.00003f;

    // The signal needs to be low for at least that many frames
    // to be uncompressed
    static constexpr int safeDuration = 100000;

    Compressor() {
#if IMJ_DEBUG_COMPRESSOR
        std::cout << "Compressor targetMultiplicator = " << targetMultiplicator
                  << " (initial value)" << std::endl;
        std::cout << "Compressor multiplicator = " << multiplicator
                  << " (initial value)" << std::endl;
#endif
    }

    template<typename T>
    void feed(T & signal) {
      Assert(targetMultiplicator <= 1.f);
      Assert(targetMultiplicator > 0.f);

      using value_type = typename T::value_type;

      value_type rawAmplitude{};
      for(auto s : signal) {
        rawAmplitude = std::max(rawAmplitude, std::abs(s));
      }

      // Adjust the _target_ compression level,
      // taking only into account the amplitude of the signal "when the compression level will have reached its target",
      if(rawAmplitude * targetMultiplicator > high) {
        // the signal is too loud
        do {
            auto targetAmplitude = rawAmplitude * targetMultiplicator;
            targetMultiplicator *= compress;
#if IMJ_DEBUG_COMPRESSOR
            std::cout << "Compressor targetMultiplicator = " << targetMultiplicator
                      << " (signal too high : targetAmplitude = " << targetAmplitude << " > " << high << " )" << std::endl;
#endif
        } while (rawAmplitude * targetMultiplicator > high);

        safeSince.onUnsafe();
      }
      else if(rawAmplitude * targetMultiplicator < low) {
        // the signal is low
        safeSince.onSafe(rawAmplitude);
        if(safeSince.getSafeSince() > safeDuration) {
          // The signal has been too low for a long time
          if (targetMultiplicator < 1.f) {
            float const maxRawAmplitude = safeSince.getSafeMaxAmplitude();
            // we adjust the target such that the maxRawAmplitude would correspond to medium level:
            //   targetMultiplicator * maxRawAmplitude = medium
            //   targetMultiplicator                   = medium / maxRawAmplitude
            if (maxRawAmplitude) {
              targetMultiplicator = std::min(1.f,medium / maxRawAmplitude);
            } else {
              targetMultiplicator = 1.f;
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
      value_type amplitude = rawAmplitude * multiplicator;
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
    float getCompressionLevel() const { return multiplicator; }
    float getTargetCompressionLevel() const { return targetMultiplicator; }
    int getSafeSince() const { return safeSince.getSafeSince(); }
  private:
    float multiplicator = 1.f;
    float targetMultiplicator = 1.f;
    SafeSince safeSince;
    const float medium = std::exp( (std::log(high) + std::log(low)) / 2.f);
  };
}
