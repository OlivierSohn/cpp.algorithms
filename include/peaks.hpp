/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule
{
    template<typename ITER, typename VAL = typename ITER::value_type>
    ITER first_relevant_value(ITER it, ITER end, VAL abs_relevant_level) {
        for(; it != end; ++ it ) {
            if(std::abs(*it) > abs_relevant_level) {
                break;
            }
        }
        return it;
    }

    template<typename ITER>
    ITER first_zero_crossing(ITER it, ITER end) {
        using VAL = typename ITER::value_type;
        using Tr = NumTraits<VAL>;

        bool first = true;
        VAL prev;
        while(it != end) {
            auto cur = *it;
            if(cur == Tr::zero()) {
                break;
            }
            if(first) {
                first = false;
            }
            else if(prev * cur < Tr::zero()) {
                break;
            }
            prev = cur;
            ++it;
        }

        return it;
    }

    template<typename ITER>
    ITER first_non_abs_decreasing(ITER it, ITER end) {
        using VAL = typename ITER::value_type;

        bool first = true;
        VAL prev;
        auto prev_it = it;
        while(it != end) {
            auto cur = std::abs(*it);
            if(first) {
                first = false;
            }
            else if(prev < cur ) {
                break;
            }
            prev = cur;
            prev_it = it;
            ++it;
        }

        return prev_it;
    }

    template<typename ITER>
    ITER first_non_abs_avg_decreasing(ITER it, ITER end, int avg_size) {
        using VAL = typename ITER::value_type;

        slidingAverage<VAL> avg(avg_size);

        bool first = true;
        VAL prev_avg;
        auto prev_it = it;
        while(it != end) {
            avg.feed(std::abs(*it));
            auto cur_avg = avg.compute();
            if(first) {
                first = false;
            }
            else if(prev_avg < cur_avg ) {
                break;
            }
            prev_avg = cur_avg;
            prev_it = it;
            ++it;
        }

        return prev_it;
    }

    template<typename ITER, typename VAL = typename ITER::value_type>
    ITER find_relevant_start(ITER it, ITER end, VAL abs_relevant_level) {
        auto it_relevant_value = first_relevant_value(it, end, abs_relevant_level);
        if(it_relevant_value == end) {
            return end;
        }
        using REVERSE_ITER = std::reverse_iterator<ITER>;
        auto rit = REVERSE_ITER(it_relevant_value + 1);
        auto rend = REVERSE_ITER(it);
        auto rzero = first_zero_crossing( rit, rend);
        // first_zero_crossing returns the iterator after the zero crossing (in the reverse direction)
        // so rzero.base() is the iterator on the other side of the zero crossing
        return rzero.base();
    }

  // returns the iterator that comes BEFORE the first relevant value.
    template<typename ITER, typename VAL = typename ITER::value_type>
    ITER find_relevant_start_relaxed(ITER const it, ITER const end, VAL const abs_relevant_level, int const sliding_avg_size) {
        auto it_relevant_value = first_relevant_value(it, end, abs_relevant_level);
        if(it_relevant_value == it) {
            return it;
        }
        if(it_relevant_value == end) {
          // not found
          return end;
        }
        using REVERSE_ITER = std::reverse_iterator<ITER>;
        auto rit = REVERSE_ITER(it_relevant_value+1);
        auto rend = REVERSE_ITER(it);
        auto rzero = first_non_abs_avg_decreasing( rit, rend, sliding_avg_size);
        // first_non_abs_avg_decreasing returns the iterator after the zero crossing (in the reverse direction)
        // so rzero.base() is the iterator on the other side of the zero crossing
        return rzero.base();
    }

    template<typename ITER, typename VAL = typename ITER::value_type>
    VAL abs_integrated(ITER it, ITER end) {
        VAL ret {};

        for(; it != end; ++it) {
            ret += std::abs(*it);
        }

        return ret;
    }

    struct FreqAmplitude {
        float relative_freq;
        float amplitude;
    };

    template<typename ITER, typename VAL = typename ITER::value_type>
    FreqAmplitude max_freq_amplitude(ITER it, ITER end) {
        using namespace fft;
        using Tag = Fastest;
        using Algo = Algo_<Tag, VAL>;
        using RealFBins = RealFBins_<Tag, VAL>;
        using ScopedContext = ScopedContext_<Tag, VAL>;

        a64::vector<VAL> v;
        auto fft_length = ceil_power_of_two(std::distance(it, end));
        v.reserve(fft_length);
        for(; it!=end; ++it) {
            v.push_back(*it);
        }
        v.resize(fft_length, VAL{0});
        auto signal = RealSignal_<Tag, VAL>::make(std::move(v));

        typename RealFBins::type result(fft_length);

        ScopedContext scoped_context(fft_length);
        Algo fft(scoped_context.get());
        fft.forward(signal, result, fft_length);
        auto Max = RealFBins::getMaxSquaredAmplitude(result);
        return {
            Max.first/static_cast<float>(fft_length),
            std::sqrt(Max.second)
        };
    }
}
