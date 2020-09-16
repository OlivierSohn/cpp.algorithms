

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


void deduceNotes() {
  int const window_center_stride = 400;
  int const windowed_signal_stride = 1;
  std::string prefix("chromatic");

  // The frequency detection is too imprecise with hann window:

  //std::vector<double> half_window = half_hann_window<double>(400);
  //std::vector<double> half_window = half_hann_window<double>(4000);
  
  // The frequency detection is very precise with a gaussian window truncated at 8 sigma
  // however, low frequencies are incorrectly detected with such a narrow window, so
  // we can either:
  // - use a window with more points (but lose some temporal precision)
  // - or truncate the window at 4 sigma for example (but lose some frequency precision)
  //
  // TODO we could use a 8-sigma window for high frequencies, and a 4-sigma window for low frequencies,
  // this way we retain temporal precision in the high frequencies, and trade temporal precision
  // for frequency accuracy only for lower frequencies, where this is needed.
  
  std::vector<double> half_window = half_gaussian_window<double>(8, 400);
  //std::vector<double> half_window = half_gaussian_window<double>(8, 4000);
  int const zero_padding_factor = 1;
  //int const zero_padding_factor = 10;
  auto reader = WAVReader("/Users/Olivier/", prefix + ".wav");

  std::vector<std::vector<double>> deinterlaced;
  read_wav_as_floats(reader, deinterlaced);

  auto notes = deduceNotes(deinterlaced.begin()->begin(),
                           deinterlaced.begin()->end(),
                           reader.getSampleRate(),
                           windowed_signal_stride,
                           half_window,
                           window_center_stride,
                           zero_padding_factor,
                           0.05776226504);
  // formula in http://support.ircam.fr/docs/AudioSculpt/3.0/co/Window%20Size.html
  // uses 5 as constant factor
  const double lowest_detectable_frequency = 4. * reader.getSampleRate() / (2*half_window.size() * windowed_signal_stride);
  std::cout << "lowest detectable freq : " << lowest_detectable_frequency << " Hz" << std::endl;

  drawDeducedNotes(notes,
                   lowest_detectable_frequency,
                   "/Users/Olivier/" + prefix + ".notes.bmp");
}

} // NS

TEST(FreqTime, drawSpectrum) {
  imajuscule::audio::drawSpectrum();
}

TEST(FreqTime, deduceNotes) {
  imajuscule::audio::deduceNotes();
}
