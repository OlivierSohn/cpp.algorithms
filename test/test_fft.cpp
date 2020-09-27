
/*
 * Comparative analysis of space complexity, for forward fft of real input of size N.
 *
 *   Summary : AccelerateFFT has 2x better cache locality
 *
 * - MyFFT :
 *
 * roots :   N
 * input : 2*N
 * result: 2*N
 *
 * - AccelerateFFT:
 *
 * setup :   N (?)
 * input :   N
 * result:   N
 *
 */

namespace imajuscule {
    namespace testfft {
        
        template<typename T>
        void createRealInput(T & input) {
            for(int i=0; i<4; ++i) {
                input.emplace_back(1);
            }
            for(int i=0; i<4; ++i) {
                input.emplace_back(0);
            }
        }
        
        template<typename T>
        auto makeCoefficients() {
            return a64::vector<T>{{ .9,.8,.7,.6,.3,.2,.1,0. }};
        }
        
        template<typename Tag, typename T, typename Input>
        void verifySignal2(Input signal, int N) {
            using namespace imajuscule::fft;
            
            const auto ffteps = getFFTEpsilon<T>(N);
            
            // we use ffteps when we have an exact value to compare with:
            
            auto const res = slow_debug::unwrap_signal<Tag>(signal, signal.size());
            
            constexpr auto scale = 1;
            int i=0;
            for(auto c : makeCoefficients<T>()) {
                ASSERT_NEAR(scale * c, res[i].real(), ffteps);
                ASSERT_NEAR(scale * 0, res[i].imag(), ffteps);
                ++i;
            }
        }

        template<typename Tag, typename T, typename Output>
        void verifyFrequencies(Output output, int N) {
            using namespace imajuscule::fft;
            
            // testing against ref found here : https://rosettacode.org/wiki/Fast_Fourier_transform
            
            const auto ffteps = getFFTEpsilon<T>(N);
            
            // we use ffteps when we have an exact value to compare with:
            
            auto const res = slow_debug::unwrap_frequencies<Tag>(output, N);
            
            constexpr auto scale = Algo_<Tag, T>::scale;
            
            ASSERT_NEAR(scale * 0, res[0].imag(), ffteps);
            ASSERT_NEAR(scale * 4, res[0].real(), ffteps);
            ASSERT_NEAR(scale * 0, res[2].imag(), ffteps);
            ASSERT_NEAR(scale * 0, res[2].real(), ffteps);
            ASSERT_NEAR(scale * 0, res[4].imag(), ffteps);
            ASSERT_NEAR(scale * 0, res[4].real(), ffteps);
            ASSERT_NEAR(scale * 0, res[6].imag(), ffteps);
            ASSERT_NEAR(scale * 0, res[6].real(), ffteps);
            
            ASSERT_NEAR(scale * 1, res[1].real(), ffteps);
            ASSERT_NEAR(scale * -2.41421, res[1].imag(), 1e-4);
            ASSERT_NEAR(scale * 1, res[3].real(), ffteps);
            ASSERT_NEAR(scale * -0.414214, res[3].imag(), 1e-4);
            ASSERT_NEAR(scale * 1, res[5].real(), ffteps);
            ASSERT_NEAR(scale * 0.414214, res[5].imag(), 1e-4);
            ASSERT_NEAR(scale * 1, res[7].real(), ffteps);
            ASSERT_NEAR(scale * 2.41421, res[7].imag(), 1e-4);
        }
        
        template<typename Tag, typename T, typename Input>
        void verifySignal(Input signal, int N) {
            using namespace imajuscule::fft;
            
            const auto ffteps = getFFTEpsilon<T>(N);
            
            // we use ffteps when we have an exact value to compare with:
                        
            auto const res = slow_debug::unwrap_signal<Tag>(signal, N);
            
            constexpr T scale = 1;
            
            ASSERT_NEAR(scale * 1, res[0].real(), ffteps);
            ASSERT_NEAR(scale * 0, res[0].imag(), ffteps);
            ASSERT_NEAR(scale * 1, res[1].real(), ffteps);
            ASSERT_NEAR(scale * 0, res[1].imag(), ffteps);
            ASSERT_NEAR(scale * 1, res[2].real(), ffteps);
            ASSERT_NEAR(scale * 0, res[2].imag(), ffteps);
            ASSERT_NEAR(scale * 1, res[3].real(), ffteps);
            ASSERT_NEAR(scale * 0, res[3].imag(), ffteps);
            ASSERT_NEAR(scale * 0, res[4].real(), ffteps);
            ASSERT_NEAR(scale * 0, res[4].imag(), ffteps);
            ASSERT_NEAR(scale * 0, res[5].real(), ffteps);
            ASSERT_NEAR(scale * 0, res[5].imag(), ffteps);
            ASSERT_NEAR(scale * 0, res[6].real(), ffteps);
            ASSERT_NEAR(scale * 0, res[6].imag(), ffteps);
            ASSERT_NEAR(scale * 0, res[7].real(), ffteps);
            ASSERT_NEAR(scale * 0, res[7].imag(), ffteps);
        }

        template<typename Tag, typename T, template<typename> typename Allocator>
        void testForwardFFT() {
            using namespace imajuscule;
            using namespace imajuscule::fft;
            using namespace imajuscule::testfft;

            using RealInput = typename RealSignal_<Tag, T>::type;
            using RealFBins = typename RealFBins_<Tag, T, Allocator>::type;
            using Context   = typename Context_<Tag, T>::type;
            using ScopedContext = ScopedContext_<Tag, T>;
            using Algo = Algo_<Tag, T>;

            constexpr auto N = 8;
            ScopedContext setup(N);
            
            RealInput input;
            RealFBins output(N);
            
            input.reserve(N);
            
            createRealInput(input);

            { verifySignal<Tag, T>(input, N); }

            Algo fft_algo(setup.get());
            
            fft_algo.forward(input.begin(), output.data(), N);
            
            { verifyFrequencies<Tag, T>(output, N); }
        }
        
        template<typename Tag, template<typename> typename Allocator>
        void testForwardFFT() {
            testForwardFFT<Tag, float, Allocator>();
            testForwardFFT<Tag, double, Allocator>();
        }
        
        template<typename Tag, typename T, template<typename> typename Allocator>
        void testInverseFFT() {
            using namespace imajuscule;
            using namespace imajuscule::fft;
            using namespace imajuscule::testfft;
            
            using RealInput = typename RealSignal_<Tag, T>::type;
            using RealFBins = typename RealFBins_<Tag, T, Allocator>::type;
            using Context   = typename Context_<Tag, T>::type;
            using ScopedContext = ScopedContext_<Tag, T>;
            using Algo = Algo_<Tag, T>;
            
            constexpr auto N = 8;
            ScopedContext setup(N);
            
            RealInput input, reconstructed_input(N);
            RealFBins output(N);
            
            input.reserve(N);
            
            createRealInput(input);
            
            { verifySignal<Tag, T>(input, N); }
            
            Algo fft_algo(setup.get());
            
            fft_algo.forward(input.begin(), output.data(), N);
            
            { verifyFrequencies<Tag, T>(output, N); }
            
            if constexpr (Algo::inplace_dft) {
                fft_algo.inverse(output.data(),
                                 N);
                for(int i=0; i<N; ++i) {
                    reconstructed_input[i] = Algo::extractRealOutput(output.data(),
                                                                     i,
                                                                     N);
                }
            }
            else {
                fft_algo.inverse(output.data(),
                                 reconstructed_input.data(),
                                 N);
            }
            
            for(auto & v : reconstructed_input) {
                v *= 1/(Algo::scale * static_cast<T>(N));
            }
            
            { verifySignal<Tag, T>(reconstructed_input, N); }
            
        }
        
        template<typename Tag, typename T, template<typename> typename Allocator>
        void testInverseFFT2() {
            using namespace imajuscule;
            using namespace imajuscule::fft;
            using namespace imajuscule::testfft;
            
            using RealInputT = RealSignal_<Tag, T>;
            using RealInput = typename RealSignal_<Tag, T>::type;
            using RealFBins = typename RealFBins_<Tag, T, Allocator>::type;
            using Context   = typename Context_<Tag, T>::type;
            using ScopedContext = ScopedContext_<Tag, T>;
            using Algo = Algo_<Tag, T>;
            
            constexpr auto N = 8;
            ScopedContext setup(N);
            
            RealInput input, reconstructed_input(N);
            RealFBins output(N);
            
            input = RealInputT::make(makeCoefficients<T>());
            EXPECT_EQ(N, input.size());
            
            { verifySignal2<Tag, T>(input, N); }
            
            Algo fft_algo(setup.get());
            
            fft_algo.forward(input.begin(), output.data(), N);
            
            if constexpr (Algo::inplace_dft) {
                fft_algo.inverse(output.data(),
                                 N);
                for(int i=0; i<N; ++i) {
                    reconstructed_input[i] = Algo::extractRealOutput(output.data(),
                                                                     i,
                                                                     N);
                }
            }
            else {
                fft_algo.inverse(output.data(),
                                 reconstructed_input.data(),
                                 N);
            }

            for(auto & v : reconstructed_input) {
                v *= 1/(Algo::scale * static_cast<T>(N));
            }
            
            { verifySignal2<Tag, T>(reconstructed_input, N); }
            
        }
        
        template<typename Tag, template<typename> typename Allocator>
        void testInverseFFT() {
            testInverseFFT<Tag, float, Allocator>();
            testInverseFFT<Tag, double, Allocator>();
            testInverseFFT2<Tag, float, Allocator>();
            testInverseFFT2<Tag, double, Allocator>();
        }
    }
}

TEST(FFT, forward_correctness) {
    using namespace imajuscule;
    using namespace imajuscule;
    using namespace imajuscule::testfft;
    
    for_each(fft::Tags, [](auto t) {
        testForwardFFT<decltype(t), a64::Alloc>();
    });
}

TEST(FFT, inverse_correctness) {
    using namespace imajuscule;
    using namespace imajuscule::testfft;
    
    for_each(fft::Tags, [](auto t) {
        testInverseFFT<decltype(t), a64::Alloc>();
    });
}

TEST(FFT, unwrap_freq_sqmag) {
  using namespace imajuscule;
  using namespace imajuscule::fft;
  using namespace imajuscule::audio;

  using Tag = fft::Fastest;
  using T = double;
  
  using RealInput = typename RealSignal_<Tag, T>::type;
  using RealFBins = typename RealFBins_<Tag, T, a64::Alloc>::type;
  using Context   = typename Context_<Tag, T>::type;
  using ScopedContext = ScopedContext_<Tag, T>;
  using Algo = Algo_<Tag, T>;
  
  constexpr int N = 1024;

  //for (int i = 0; i<=10; ++i)
  int const i = 10;
  for (int M = 1; M <= 1024; M/**=2*/++)
  // for M=1 and 2, the signal is almost always zero.
  // for low frequencies (M between 513 and 1024), the amplitude detection error is huge
  {
    ScopedContext setup(N);
    
    RealInput input0, input;
    RealFBins output(N);
    
    input.reserve(N);
    
    double const amplitude = i / 10.f;
    std::cout << "amplitude " << amplitude << " M=" << M << std::endl;
    for (int i=0; i<N; ++i) {
      input.push_back(amplitude * std::sin(2. * M_PI * static_cast<double>(i) / M));
    }
    
    std::vector<double> half_window = half_gaussian_window<double>(4, N/2);
    //for (auto & v : half_window) {v = 1.;}
    normalize_window(half_window);
    
    int const signal_stride = 1;
    apply_window(input.begin(), input.end(), signal_stride, half_window, input);
    
    Algo fft_algo(setup.get());
    
    fft_algo.forward(input.begin(), output.data(), N);
    
    FrequenciesSqMag<double> frequencies_sqmag;
    unwrap_frequencies_sqmag<Tag>(output, N, frequencies_sqmag.frequencies_sqmag);
    frequencies_sqmag.fft_length = N;

    constexpr T inv_scale_squared = 1. / (Algo::scale * Algo::scale * N * N);
    for (auto & v : frequencies_sqmag.frequencies_sqmag) {
      v *= inv_scale_squared;
    }
    
    double const sample_rate = N;
    std::vector<FreqMag<double>> freqmags;
    extractLocalMaxFreqsMags(sample_rate / signal_stride,
                             frequencies_sqmag,
                             SqMagToDb<double>(),
                             freqmags);
    
    for (auto const & r : freqmags) {
      std::cout << r << " mag= " << DbToMag<double>()(r.mag_db) << std::endl;
    }

    std::optional<std::pair<int, T>> maxValue;
    int idx = -1;
    for (auto const & m : frequencies_sqmag.frequencies_sqmag) {
      ++idx;
      if (!maxValue || maxValue->second < m) {
        maxValue = std::make_pair(idx, m);
      }
    }
    
    for (int i = std::max(0, maxValue->first-2);
         i<=std::min(maxValue->first+2,
                     static_cast<int>(frequencies_sqmag.frequencies_sqmag.size())-1);
         ++i) {
      std::cout << i << ":" << std::sqrt(frequencies_sqmag.frequencies_sqmag[i]) << std::endl;
    }
  }
}
