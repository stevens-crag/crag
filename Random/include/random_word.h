#pragma once

#ifndef CRAG_RANDOM_WORD_H
#define CRAG_RANDOM_WORD_H

#include <algorithm>
#include <boost/random/uniform_int_distribution.hpp>

#include "Permutation.h"
#include "Word.h"

namespace crag {
namespace random {

//! Generates pseudo-randomly a reduced word of the length len over the alphabet \f$\{x_1, \ldots, x_n\}\f$.
template <typename URNG>
Word randomWord(size_t n, size_t len, URNG& g) {
  if (len == 0) {
    return Word();
  }

  std::list<int> result;

  const auto dist_first = boost::random::uniform_int_distribution<int>(-n, n - 1);
  const auto dist_rest = boost::random::uniform_int_distribution<int>(-n, n - 2);

  auto random_number = dist_first(g);

  if (random_number >= 0) {
    random_number += 1;
  }

  result.push_back(random_number);

  auto old_random_number = random_number;

  for (size_t i = 1; i < len; ++i) {
    random_number = dist_rest(g);

    if (random_number >= 0) {
      random_number += 1;
    }

    if (random_number + old_random_number == 0) {
      random_number = n;
    }

    result.push_back(random_number);

    old_random_number = random_number;
  }

  return Word(std::move(result));
}

//! Generates pseudo-randomly a reduced word of a length in [from, to] over the alphabet \f$\{x_1, \ldots, x_n\}\f$.
template <typename URNG>
Word randomWord(size_t n, size_t from, size_t to, URNG& g) {
  const auto len = boost::random::uniform_int_distribution<>(from, to)(g);
  return randomWord(n, len, g);
}

//! Generates a random braid word with a given permutation.
//! Length of the word is determined by the length of p in (i, i+1) generators.
template <typename URNG>
Word randomWord(const Permutation& p, URNG& g) {
  // First write p as a product of disjoint cycles C_1,...,C_s
  // where the last element of each C_i is the smallest number in the cycle
  auto cycles = toCycles(p);

  for (auto& c : cycles) {
    auto it = std::min_element(c.begin(), c.end());

    if (++it != c.end()) {
      std::rotate(c.begin(), it, c.end());
    }
  }

  // Order the cycles such that the last element of each C_i is in ascending order
  std::sort(cycles.begin(), cycles.end(), [](const std::vector<int>& lhs, const std::vector<int>& rhs) {
    return lhs.back() < rhs.back();
  });

  Word w;

  boost::random::uniform_int_distribution<int> dist(0, 1);

  for (const auto& c : cycles) {
    for (size_t i = 1; i < c.size(); ++i) {
      const auto m = std::min(c[0], c[i]);
      const auto M = std::max(c[0], c[i]);

      // generate word for transposition (m, M)
      for (int k = m; k < M; ++k) {
        // append b_{k + 1}^{\pm 1}
        w.push_back((1 - 2 * dist(g)) * (k + 1));
      }

      for (int k = 2; k <= M - m; ++k) {
        // append b_{M - k + 1}^{\pm 1}
        w.push_back((1 - 2 * dist(g)) * (M - k + 1));
      }
    }
  }

  return w;
}
} // namespace random
} // namespace crag

#endif /* CRAG_RANDOM_WORD_H */
