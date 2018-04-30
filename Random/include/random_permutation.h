#pragma once

#ifndef CRAG_RANDOM_PERMUTATION_H
#define CRAG_RANDOM_PERMUTATION_H

#include <numeric>

#include "Permutation.h"
#include "shuffle.h"

namespace crag {
namespace random {

template <typename URNG>
Permutation randomPermutation(size_t n, URNG& g) {
  std::vector<int> values(n);
  std::iota(values.begin(), values.end(), 0);

  random::shuffle(values.begin(), values.end(), g);

  return Permutation(std::move(values));
}
} // namespace random
} // namespace crag

#endif // CRAG_RANDOM_PERMUTATION_H
