/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule
{
// Returns minimum of "max abolute value" over a range of 'interval_size' consecutive values.
template<typename ITER, typename VAL = typename ITER::value_type>
VAL compute_noise_floor(ITER it, ITER end, int interval_size)
{  
  slidingAverage<VAL> avg(interval_size);

  VAL max{std::numeric_limits<VAL>::max()};

  int i = 0;
  for(;it != end; ++it, ++i) {
    avg.feed(std::abs(*it));
    if(i >= interval_size)
    {
      auto cur_max = avg.computeMax();
      max = std::min(max, cur_max);
    }
  }
  
  return max;
}

// Returns the first iterator that has an absolute value larger than |abs_relevant_level|.
// Otherwise, returns end.
template<typename ITER, typename VAL = typename ITER::value_type>
ITER first_relevant_value(ITER it, ITER end, VAL abs_relevant_level) {
  for(; it != end; ++ it ) {
    if(std::abs(*it) > abs_relevant_level) {
      break;
    }
  }
  return it;
}

// Returns the first iterator that has a 0 value or that has a sign change wrt the previous iterator.
// Otherwise, returns end.
//
// Note: it is used with reverse iterators to find the beginning of a sample,
//       once we have found the approximate beginning of the sample using first_relevant_value.
template<typename ITER>
ITER first_zero_crossing(ITER it, ITER end) {
  using VAL = typename ITER::value_type;
  using Tr = NumTraits<VAL>;
  
  std::optional<VAL> prev;
  while(it != end) {
    auto cur = *it;
    if(cur == Tr::zero()) {
      break;
    }
    if(prev && (*prev * cur < Tr::zero())) {
      break;
    }
    prev = cur;
    ++it;
  }
  
  return it;
}

// Returns the first iterator that has an absolute value that increases wrt the previous iterator.
// Otherwise returns end.
template<typename ITER>
ITER first_non_abs_decreasing(ITER it, ITER end) {
  using VAL = typename ITER::value_type;
  
  std::optional<VAL> prev;
  auto prev_it = it;
  while(it != end) {
    auto cur = std::abs(*it);
    if(prev && (*prev < cur)) {
      break;
    }
    prev = cur;
    prev_it = it;
    ++it;
  }
  
  return prev_it;
}

// Using a sliding average of absolute values, returns the first iterator
// for which the sliding average increases wrt the previous iterator.
// Otherwise returns end.
//
// Typical sliding average length : 15
template<typename ITER>
ITER first_non_abs_avg_decreasing(ITER it, ITER end, int avg_size) {
  using VAL = typename ITER::value_type;
  
  slidingAverage<VAL> avg(avg_size);
  
  std::optional<VAL> prev_avg;
  auto prev_it = it;
  while(it != end) {
    avg.feed(std::abs(*it));
    auto cur_avg = avg.compute();
    if(prev_avg && (*prev_avg < cur_avg)) {
      break;
    }
    prev_avg = cur_avg;
    prev_it = it;
    ++it;
  }
  
  return prev_it;
}

// Moves fwd until the sliding average of absolute signal becomes lower than
//     'abs_relevant_level' for at least 'lookahead_size' frames.
// Returns the location where the signal started to be low.
template<typename ITER, typename VAL = typename ITER::value_type>
ITER first_avg_non_relevant_value(ITER it, ITER end,
                                  VAL const abs_relevant_level,
                                  int avg_size,
                                  int32_t const lookahead_size)
{
  slidingAverage<VAL> avg(avg_size);

  ITER res = end;

  int32_t countMatches{};
  for(;it != end; ++it) {
    avg.feed(std::abs(*it));
    auto cur_avg = avg.compute();
    if(cur_avg <= abs_relevant_level) {
      if(countMatches == 0)
        res = it;
      ++countMatches;
      if(countMatches > lookahead_size)
        break;
    }
    else
      countMatches = 0;
  }
  if(countMatches)
    return res;
  return end;
}

// Uses first_relevant_value(fwd, abs_relevant_level) then first_zero_crossing(bwd)
//
// returns the iterator that comes BEFORE the first relevant value.
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

// Uses first_relevant_value(fwd, abs_relevant_level) then first_non_abs_avg_decreasing(bwd, sliding_avg_size)
//
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

// find the end of a range:
// - move fwd until signal (sliding average of absolute) is lower than
//     'abs_relevant_level' for at least 'lookahead_size' frames.
// - from first location where signal started to be low,
//   move fwd until (sliding average of absolute) increases again.
//
// returns the iterator that is the last relevant value,
// or end if all values are relevant.
template<typename ITER, typename VAL = typename ITER::value_type>
ITER find_relevant_end_relaxed(ITER const it, ITER const end,
                               VAL const abs_relevant_level,
                               int const sliding_avg_size,
                               int32_t const lookahead_size) {
  auto it_non_relevant_value = first_avg_non_relevant_value(it, end, abs_relevant_level, sliding_avg_size, lookahead_size);
  if(it_non_relevant_value == end) {
    return end;
  }
  return first_non_abs_avg_decreasing(it_non_relevant_value, end, sliding_avg_size);
}


template<typename ITER, typename VAL = typename ITER::value_type>
VAL abs_integrated(ITER it, ITER end) {
  VAL ret {};
  
  for(; it != end; ++it) {
    ret += std::abs(*it);
  }
  
  return ret;
}
}
