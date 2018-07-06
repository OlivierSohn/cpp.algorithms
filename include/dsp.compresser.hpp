
namespace imajuscule::audio {

  /* This compressor quickly reduces the signal volume
   * when the signal becomes too loud. and slowly restores
   * the initial volume when the signal has been below a threshold for a
   * significant period of time.
   *
   * To avoid having having compression level changes eventhough the signal
   * didn't change much, we use hysteresis.
   */
  struct Compressor {
    // compress == 0.5f was making a too big change in volume.
    static constexpr float compress   = 0.8f;
    static constexpr float uncompress = 1/compress;

    static constexpr float high = 0.7f;
    static constexpr float low  = 0.6f * compress;

    static constexpr float compressMult   = 0.997f;

    // = 2^(1/22050) so that the volume doubles every half second
    static constexpr float uncompressMult = 1.00003f;


    static constexpr int safeDuration = 100000;

    template<typename T>
    void feed(T & signal) {
      using value_type = typename T::value_type;
      // find the max amplitude
      value_type targetAmplitude{};
      for(auto preS : signal) {
        // s is the value of the signal when the compression level will have reached its target.
        auto s = preS * targetMultiplicator;
        targetAmplitude = std::max(targetAmplitude, std::abs(s));
      }

      // adjust the target compression level
      if(targetAmplitude > high) {
        // the signal is too loud
        targetMultiplicator *= compress;
        safeSince = 0;
      }
      else if(targetAmplitude < low) {
        // the signal is low
        ++safeSince;
        if(safeSince > safeDuration) {
          targetMultiplicator = std::min(1.f,targetMultiplicator*uncompress);
          safeSince = 0;
        }
      }
      else {
        // the signal is normal, compression is good.
        safeSince = 0;
      }

      // adjust the compression level
      if(likely(multiplicator == targetMultiplicator)) {
        // do nothing
      }
      else {
        if(multiplicator < targetMultiplicator) {
          multiplicator = std::min(multiplicator * uncompressMult, targetMultiplicator);
        }
        else {
          multiplicator = std::max(multiplicator * compressMult, targetMultiplicator);
        }
      }

      // compress the signal.
      for(auto & s : signal) {
        s *= multiplicator;
      }
    }

    // needed for white box tests
    float getCompressionLevel() const { return multiplicator; }
    float getTargetCompressionLevel() const { return targetMultiplicator; }
    int getSafeSince() const { return safeSince; }
  private:
    float multiplicator = 1.f;
    float targetMultiplicator = 1.f;
    // the multiplicator can go up, down, or be stationnary.
    int safeSince = 0;

  };
}
