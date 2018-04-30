#pragma once

#ifndef CRAG_RANDOM_SUBSET_H
#define CRAG_RANDOM_SUBSET_H

#include <boost/random/uniform_int_distribution.hpp>
#include <set>
#include <type_traits>
#include <vector>

#include "shuffle.h"

namespace crag {
namespace random {

//! Returns a random n-element subset of integers from [min, max].
template <typename T, typename URNG, typename = typename std::enable_if<std::is_integral<T>::value>::type>
std::set<T> subset(size_t n, T min, T max, URNG& g) {
  if (max < min) {
    throw std::invalid_argument("Maximal value must be greater than or equal to minimal.");
  }

  if ((max - min) + 1 < n) {
    throw std::invalid_argument("The size of the desired random subset is too large.");
  }

  boost::random::uniform_int_distribution<T> distribution(min, max);

  std::set<T> result;

  while (result.size() < n) {
    const auto value = distribution(g);

    if (result.count(value) == 0) {
      result.insert(value);
    }
  }

  return result;
}

//! Returns a random n-element subset of integers from [min, max] as an unordered randomly shuffled vector.
template <typename T, typename URNG, typename = typename std::enable_if<std::is_integral<T>::value>::type>
std::vector<T> subsetVector(size_t n, T min, T max, URNG& g) {
  const auto s = subset(n, min, max, g);

  std::vector<T> result(s.begin(), s.end());

  random::shuffle(result.begin(), result.end(), g);

  return result;
}
} // namespace random
} // namespace crag

#endif // CRAG_RANDOM_SUBSET_H
