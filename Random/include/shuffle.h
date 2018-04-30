#pragma once

#ifndef CRAG_SHUFFLE_H
#define CRAG_SHUFFLE_H

#include <boost/random/uniform_int_distribution.hpp>

namespace crag {
namespace random {

//! Shuffles elements in a given range,
//! Workaround for std::shuffle which uses std::uniform_int_distribution,
//! which is inconsistent across different platforms.
template <typename RandomAccessIterator, typename URNG>
void shuffle(RandomAccessIterator first, RandomAccessIterator last, URNG&& g) {
  using diff_t = typename std::iterator_traits<RandomAccessIterator>::difference_type;
  using distr_t = boost::random::uniform_int_distribution<diff_t>;
  using param_t = typename distr_t::param_type;

  distr_t distribution;
  diff_t n = last - first;

  for (diff_t i = n - 1; i > 0; --i) {
    std::swap(first[i], first[distribution(g, param_t(0, i))]);
  }
}
} // namespace random
} // namespace crag

#endif // CRAG_SHUFFLE_H
