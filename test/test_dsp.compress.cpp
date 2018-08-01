
TEST(DspCompress, simple) {
  using namespace imajuscule;
  using namespace imajuscule::audio;
  Compressor c;
  
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
    ASSERT_FLOAT_EQ(Compressor::compress, c.getTargetCompressionLevel());
  }
  ASSERT_FLOAT_EQ(Compressor::compress, c.getCompressionLevel());
  
  // a 0.69 signal, just below the upper limit, should not trigger uncompression.
  for(int i=0; i<100000; ++i) {
    std::array<float,1> signal{{0.69f}};
    c.feed(signal);
    ASSERT_EQ(0, c.getSafeSince());
    ASSERT_FLOAT_EQ(Compressor::compress, c.getTargetCompressionLevel());
    ASSERT_FLOAT_EQ(Compressor::compress, c.getCompressionLevel());
    
    ASSERT_FLOAT_EQ(Compressor::compress * 0.69f, signal[0]);
  }
  
  // a 0.59 signal, just below the lower limit, should trigger uncompression (after some time).
  for(int i=0; i<Compressor::safeDuration; ++i) {
    std::array<float,1> signal{{0.59f}};
    c.feed(signal);
    ASSERT_EQ(i+1, c.getSafeSince());
    ASSERT_FLOAT_EQ(Compressor::compress, c.getTargetCompressionLevel());
    ASSERT_FLOAT_EQ(Compressor::compress, c.getCompressionLevel());
  }
  {
    std::array<float,1> signal{{0.59f}};
    c.feed(signal);
    ASSERT_EQ(0, c.getSafeSince());
    ASSERT_FLOAT_EQ(1.f, c.getTargetCompressionLevel());
  }
  for(int i=0; i<100000; ++i) {
    std::array<float,1> signal{{0.59f}};
    c.feed(signal);
    ASSERT_FLOAT_EQ(1.f, c.getTargetCompressionLevel());
  }
  ASSERT_FLOAT_EQ(1.f, c.getCompressionLevel());
  
}
