

namespace imajuscule::audio {


void drawSpectrumSlow() {
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
  
  drawSpectrumSlow(deinterlaced.begin()->begin(),
                   deinterlaced.begin()->end(),
                   windowed_signal_stride,
                   half_window,
                   window_center_stride,
                   zero_padding_factor,
                   "/Users/Olivier/chromatic.spectrum.bmp");
}

struct FreqComponent {
  double freq;
  double amplitude;
  double phase;
  
  double getAtTime(double t) const
  {
    return amplitude * std::sin(2. * M_PI * freq * t + phase);
  }
};

template<typename T>
void signalWithFreqAmps(std::vector<FreqComponent> const & specs,
                        int const sample_rate,
                        int const count_frames,
                        a64::vector<T> & signal) {
  signal.clear();
  signal.reserve(count_frames);
  double time_increment = 1. / sample_rate;
  for (int i=0; i < count_frames; ++i) {
    double val{};
    double time = time_increment * i;
    for (auto const & spec : specs) {
      val += spec.getAtTime(time);
    }
    signal.push_back(val);
  }
}

} // NS

TEST(FreqTime, drawSpectrumSlow) {
  imajuscule::audio::drawSpectrumSlow();
}

// This is an attempt to prototype a peak finding method where
// we iteratively substract
TEST(ParabollaGaussian, test) {
  using namespace imajuscule;
  using namespace imajuscule::audio;

  using T = double;
  using Tag = fft::Fastest;

  for (int exp = 9; exp < 20; ++exp) {
    for (int sigma = 4; sigma <= 10; ++sigma) {
      int const signal_size = pow2(exp);
      
      a64::vector<double> signal;
      typename fft::RealSignal_<Tag, T>::type work_signal;
      typename fft::RealFBins_<Tag, T, a64::Alloc>::type result;
            
      int const sample_rate = 44100;
      signalWithFreqAmps(
                         {
        {445., 0.5, 0.},
        {485., 0.5, 0.},
      },
                         sample_rate,
                         signal_size,
                         signal);
      
      // window it by a gaussian, compute spectrum, extract squared magnitudes
      
      std::vector<T> half_window;
      half_gaussian_window<T>(sigma, signal_size/2, half_window);
      normalize_window(half_window);
      
      FrequenciesSqMag<T> frequencies_sqmag;
      
      findFrequenciesSqMag<Tag>(signal.begin(),
                                signal.end(),
                                1,
                                half_window,
                                1,
                                work_signal,
                                result,
                                frequencies_sqmag);
      double const bin_index_to_Hz = frequencies_sqmag.bin_index_to_Hz(sample_rate);

      std::cout << "sigma = " << sigma << " signal size = " << signal_size << " fft size = " << frequencies_sqmag.get_fft_length() << std::endl;
            
      // verify the fit of "dB scale" amplitudes with a parabolla
      a64::vector<double> db_amplitude, amplitude;
      amplitude.clear();
      amplitude.reserve(frequencies_sqmag.frequencies_sqmag.size());
      db_amplitude.clear();
      db_amplitude.reserve(frequencies_sqmag.frequencies_sqmag.size());
      for (auto const & sqmag : frequencies_sqmag.frequencies_sqmag) {
        amplitude.push_back(std::sqrt(sqmag));
        if(sqmag <= 0) {
          db_amplitude.push_back(std::numeric_limits<double>::lowest());
        } else {
          db_amplitude.push_back(*SqMagToDb<double>()(sqmag));
        }
      }

      foreachLocalMaxFreqsMags(db_amplitude,
                               [](T value) -> std::optional<T> { return {value}; },
                               [&db_amplitude, bin_index_to_Hz](T const frequency_bin, T const magnitude){
        double const x0 = frequency_bin;
        double const y0 = magnitude;
        
        //std::cout << "peak at bin " << x0 << " bin size Hz = " << bin_index_to_Hz << std::endl;
        
        double const found_freq = bin_index_to_Hz * x0;
        double found_amplitude = DbToMag<double>()(y0);
        std::cout << "peak at freq " << found_freq << " amplitude " << found_amplitude << std::endl;
        
      });

      std::cout << "iterative substraction" << std::endl;

      // iteratively : find peak, remove peak using parabolla
      while(true) {
        std::optional<double> max;
        std::optional<int> max_index;
        for (int i=0, sz=db_amplitude.size(); i < sz; ++i) {
          if (i==0 || i==sz-1) {
            continue;
          }
          if (db_amplitude[i] == std::numeric_limits<double>::lowest()) {
            continue;
          }
          if (db_amplitude[i-1] == std::numeric_limits<double>::lowest()) {
            continue;
          }
          if (db_amplitude[i+1] == std::numeric_limits<double>::lowest()) {
            continue;
          }
          if (db_amplitude[i] >= db_amplitude[i-1] &&
              db_amplitude[i] > db_amplitude[i+1]) {
            if (!max || *max < db_amplitude[i]) {
              max = db_amplitude[i];
              max_index = i;
            }
          }
        }
        if(!max_index) {
          break;
        }

        auto itp = QuadraticInterpolation(db_amplitude[*max_index-1],
                                          db_amplitude[*max_index],
                                          db_amplitude[*max_index+1]);
        
        /* The parabolla has an equation:
         y-y0 = a*(x-x0)^2
         We know the parabolla summit (using itp) so we know 'y0' and 'x0'
         To find 'a' we use a known point which is distinct from the summit, and use:
         a = (yPoint-y0)/((xPoint-x0)^2)
         */
        double const x0 = *max_index + itp.getP();
        double const y0 = itp.getMag();
        double const xPoint = *max_index-1;
        double const yPoint = db_amplitude[*max_index-1];
        Assert(xPoint != x0);
        double const a = (yPoint-y0)/((xPoint-x0)*(xPoint-x0));

        //std::cout << "peak at bin " << x0 << " bin size = " << bin_index_to_Hz << std::endl;
        
        double const found_freq = bin_index_to_Hz * x0;
        double found_amplitude = DbToMag<double>()(y0);
        std::cout << "peak at freq " << found_freq << " amplitude " << found_amplitude << std::endl;
        
        // TODO reduce the scope to places where the parabolla is significant
        for (int i=0, sz=db_amplitude.size(); i < sz; ++i) {
          double const y = y0 + a*(i-x0)*(i-x0);
          double const new_amplitude = amplitude[i] - DbToMag<double>()(y);
          if (new_amplitude <= 0) {
            db_amplitude[i] = std::numeric_limits<double>::lowest();
          } else {
            db_amplitude[i] = *MagToDb<double>()(new_amplitude);
          }
          amplitude[i] = std::max(0., new_amplitude);
        }
      }
      //*/
    }
  }
}
