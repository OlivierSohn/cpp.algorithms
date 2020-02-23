/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule {
namespace fft {

/*
 on OSX, 'accelerate' fft is 20x faster than 'imj' fft
 */

constexpr std::tuple<
#if __APPLE__
#  ifndef IMJ_USE_SLOW_FFT
accelerate2::Tag,  // inplace
accelerate::Tag,   // osx / ios only, vectorized, fast, small memory footprint
#  endif
#endif // __APPLE__
imj::Tag,  // cross-platform, 'straightforward', slow, large memory footprint
imj2::Tag  // cross-platform, faster, smaller memory footprint
> Tags;

        using Fastest =
#ifdef IMJ_USE_SLOW_FFT
        imj::Tag;
#elif __APPLE__
        accelerate::Tag;
#else
        imj::Tag;
#endif

        template<typename T>
        constexpr double getFFTEpsilon(int N) {
            return power_of_two_exponent(N) * std::numeric_limits<T>::epsilon(); // worst case error propagation is O(log N)
        }

        template<typename T>
        using FFTVec = typename std::vector<complex<T>>;

        template<typename T>
        using FFTIter = typename FFTVec<T>::iterator;

        template<typename CONTAINER1, typename CONTAINER2>
        void forward_fft(unsigned int fft_length,
                         CONTAINER1 & signal,
                         CONTAINER2 & result) {
            using T = typename CONTAINER1::value_type;
            using FPT = typename T::FPT;

            using namespace imj::fft;

            ScopedContext<FPT> scoped_context(fft_length);
            Algo<FPT> fft(scoped_context.get());
            fft.forward(signal.begin(), result.data(), fft_length);
        }

        template<typename T>
        struct FPT {
            using type = typename T::FPT;
        };

        template<>
        struct FPT<double> {
            using type = double;
        };

        template<>
        struct FPT<float> {
            using type = float;
        };

        template<typename ITER>
        void apply_hann_window(ITER it,
                               ITER end)
        {
            using V = typename ITER::value_type;
            using T = typename FPT<V>::type;
            auto NumTaps = std::distance(it, end);
            auto nyquist = NumTaps / 2;

            int i=0;
            for(; it != end; ++it) {
                // W(n) = cos(n/NumTaps · π/2)
                *it *= std::cos( std::abs(nyquist-i)/static_cast<T>(nyquist) * static_cast<T>(M_PI_2));
                ++i;
            }
        }

    } // NS fft
} // NS imajuscule
