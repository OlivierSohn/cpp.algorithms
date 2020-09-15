

namespace imajuscule::audio {

void testMaxMagFrequencies() {
  //std::vector<double> signal {1., 1., 1., 1., -1., -1., -1., -1., 1., 1.}; // 0 (with aliazing hence sometimes 1)
  //std::vector<double> signal {1., 1., -1., -1., 1., 1., -1., -1., 1., 1.}; // 1, sq_mag 8
  std::vector<double> signal {1., -1., 1., -1., 1., -1., 1., -1., 1., -1.}; // 2, sq_mag 16
  std::vector<double> half_window{1.0, 1.0};
  std::vector<FreqSqMag> max_freqs, expected;
  for (int i=0; i<7; ++i) {
    float amp = std::log(std::sqrt(16.f));
    expected.push_back({2, amp});
  }
  findMaxMagFrequencies(signal.begin(),
                        signal.end(),
                        1,
                        half_window,
                        1,
                        max_freqs);
  EXPECT_EQ(expected, max_freqs);
}


void drawSpectrum() {
  int const window_center_stride = 400;
  int const windowed_signal_stride = 10;
  // Using a small half-window size will not be enough to detect low frequencies with a good resolution.
  // So we can use huge window sizes,

  // TODO use a multi- approach : use the same window size (eq to window_center_stride) but
  // use different signal strides (use larger signal stride + low pass) to detect lower frequencies.
  // This will be computationally less intensive (smaller ffts)

  // TODO implement low pass on signal when sigal_stride is not 1, to avoid aliasing effects

  // TODO when detecting frequency peaks, interpolate so that
  //
  //  *                        *
  //  **  is not the same as  **
  // ***                      ***
  //
  // in https://www.dsprelated.com/freebooks/sasp/Sinusoidal_Peak_Interpolation.html#:~:text=Parabolic%20interpolation%20is%20unbiased%20when,window%20transform%20about%20its%20midpoint). :
  // a ``perceptually ideal'' spectral interpolation method that is even more efficient is to zero-pad by some small factor (usually less than 5), followed by quadratic interpolation of the spectral magnitude
  //
  // and :
  // the Gaussian window transform magnitude is precisely a parabola on a dB scale. As a result, quadratic spectral peak interpolation is exact under the Gaussian window
  // ... so we should use a gaussian window for more precise quadratic quadratic spectral interpolation
  //
  // also:
  // https://ccrma.stanford.edu/~jos/parshl/Peak_Detection_Steps_3.html
  //
  // and:
  // https://www.bitweenie.com/listings/fft-zero-padding/
  //
  // TODO use gaussian window for better accuracy (https://mgasior.web.cern.ch/pap/FFT_resol_note.pdf)

  std::vector<float> half_window = half_hann_window<float>(400);
  int const zero_padding_factor = 1;
  auto reader = WAVReader("/Users/Olivier/", "melodie.wav");
  drawSpectrum(reader,
               windowed_signal_stride,
               half_window,
               window_center_stride,
               zero_padding_factor,
               "/Users/Olivier/melodie.spectrum.bmp");
}

} // NS

TEST(FreqTime, basic) {
  imajuscule::audio::testMaxMagFrequencies();
}

TEST(FreqTime, drawSpectrum) {
  imajuscule::audio::drawSpectrum();
}
