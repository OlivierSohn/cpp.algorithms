

namespace imajuscule::audio {


void drawSpectrum() {
  int const window_center_stride = 400;
  int const windowed_signal_stride = 1;

  // if we add zero padding we will have a better fft resolution

  // if we use larger windows (larger ffts) we will have a better frequency resolution

  // TODO investigate multiscale : use the same window size but on subsampled versions of the signal (resampleSinc).
  // (we wil lose some temporal precision, see if it is ok)

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
  auto reader = WAVReader("/Users/Olivier/", "chromatic.wav");
  std::vector<std::vector<float>> deinterlaced;
  read_wav_as_floats(reader, deinterlaced);

  drawSpectrum(deinterlaced.begin()->begin(),
               deinterlaced.begin()->end(),
               windowed_signal_stride,
               half_window,
               window_center_stride,
               zero_padding_factor,
               "/Users/Olivier/chromatic.spectrum.bmp");
}


} // NS

TEST(FreqTime, drawSpectrum) {
  imajuscule::audio::drawSpectrum();
}
