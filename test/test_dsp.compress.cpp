
TEST(DspCompress, simple) {
  using namespace imajuscule;
  using namespace imajuscule::audio;
  Limiter<float> c;

  ASSERT_FLOAT_EQ(1.f, c.getCompressionLevel());
  ASSERT_FLOAT_EQ(1.f, c.getTargetCompressionLevel());
  ASSERT_EQ(0, c.getSafeSince());

  // a 0 signal should not trigger compression.
  for(int i=0; i<100000; ++i) {
    std::array<float,1> signal{{0.f}};
    c.feed(signal);
    ASSERT_EQ(i + 1, c.getSafeSince());
    ASSERT_FLOAT_EQ(1.f, c.getTargetCompressionLevel());
    ASSERT_FLOAT_EQ(1.f, c.getCompressionLevel());

    ASSERT_FLOAT_EQ(0.f, signal[0]);
  }


  // a 0.69 signal, just below the upper limit, should not trigger compression.
  for(int i=0; i<100000; ++i) {
    std::array<float,1> signal{{0.69f}};
    c.feed(signal);
    ASSERT_FLOAT_EQ(1.f, c.getTargetCompressionLevel());
    ASSERT_FLOAT_EQ(1.f, c.getCompressionLevel());

    ASSERT_FLOAT_EQ(0.69f, signal[0]);
  }

  // a 0.71 signal, just above the upper limit, should trigger compression.
  for(int i=0; i<100000; ++i) {
    std::array<float,1> signal{{0.71f}};
    c.feed(signal);
    ASSERT_EQ(0, c.getSafeSince());
    ASSERT_GT(1.f, c.getTargetCompressionLevel());
  }
  ASSERT_GT(1.f, c.getCompressionLevel());

  auto const prevTarget = c.getTargetCompressionLevel();

  // a 0.69 signal, just below the upper limit, should not trigger uncompression.
  for(int i=0; i<100000; ++i) {
    std::array<float,1> signal{{0.69f}};
    c.feed(signal);
    ASSERT_EQ(0, c.getSafeSince());
    ASSERT_FLOAT_EQ(prevTarget, c.getTargetCompressionLevel());
  }

  // a 0.59 signal, just below the lower limit, should trigger uncompression (after some time)...
  for(int i=0; i<Limiter<float>::safeDuration; ++i) {
    std::array<float,1> signal{{0.59f}};
    c.feed(signal);
    ASSERT_EQ(i+1, c.getSafeSince());
    ASSERT_FLOAT_EQ(prevTarget, c.getTargetCompressionLevel());
  }
  // The current implementation has hysteresis on the target compression level
  // so 1.f will not be reached. While the signal is in the "safe" zone,
  // the compressor tracks the max amplitude, and when it is time to uncompress,
  // it will compute a target such that the compressed max amplitude will be
  // at the "medium" level (between high and low). Hence:
  for(int i=0; i<100000; ++i) {
    std::array<float,1> signal{{0.59f}};
    c.feed(signal);
    ASSERT_EQ(0, c.getSafeSince());
    ASSERT_LT(prevTarget, c.getTargetCompressionLevel());
  }
}

TEST(DspCompress, clip) {
  using namespace imajuscule::audio;
  Limiter<float> c;

  ASSERT_FLOAT_EQ(1.f, c.getCompressionLevel());
  ASSERT_FLOAT_EQ(1.f, c.getTargetCompressionLevel());
  ASSERT_EQ(0, c.getSafeSince());

  // a clipping signal (-10.) should trigger immediate compression
  std::array<float,1> signal{{-10.f}};
  c.feed(signal);
  ASSERT_EQ(0, c.getSafeSince());
  ASSERT_GT(c.getCompressionLevel(), c.getTargetCompressionLevel());
  ASSERT_FLOAT_EQ(Limiter<float>::limit / 10.f, c.getCompressionLevel());
}
