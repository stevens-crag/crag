#pragma once

#ifndef CRAG_CLOAKING_ELEMENT_H
#define CRAG_CLOAKING_ELEMENT_H

#include "braid_group.h"
#include "colored_burau.h"
#include "random_permutation.h"
#include "random_subset.h"
#include "random_word.h"

namespace crag {
namespace coloredburau {
namespace internal {

// Generate a cloaking element = a pure braid of the form w*x_j^2*w^-1, where w
// is a random braid word of length [w_min_length, w_max_length] with added suffix to satisfy
// sigma_w(a)=j-1 and sigma_w(b)=j
template <typename URNG>
Word generateCloakingElement(
    size_t n, size_t a, size_t b, size_t j, size_t w_min_length, size_t w_max_length, int center_degree, URNG& g) {
  auto w = random::randomWord(n - 1, w_min_length, w_max_length, g);
  // Make sure w satisfies cloaking condition: p[a] == j - 1 && p[b] == j
  Permutation p1 = coloredburau::permutation(n, w);
  Permutation p2 = -p1;

  while (p2[j - 1] != a || p2[j] != b) {
    while (p2[j - 1] != a) {
      const int preimage = p2[j - 1];
      int gen = (preimage < a ? preimage + 1 : preimage);
      w.push_front(g() % 2 == 0 ? gen : -gen);
      p1.change(gen - 1, gen);
      p2 = -p1;
    }

    while (p2[j] != b) {
      const int preimage = p2[j];
      int gen = (preimage < b ? preimage + 1 : preimage);
      w.push_front(g() % 2 == 0 ? gen : -gen);
      p1.change(gen - 1, gen);
      p2 = -p1;
    }
  }

  Word center = Word(j) ^ center_degree;

  if (g() % 2 != 0) {
    center = -center;
  }

  return w * center * -w;
}

template <typename URNG>
Permutation generateSigmaW(size_t n, size_t a, size_t b, size_t i, const Permutation& sigma, URNG& g) {
  auto sigma_w = random::randomPermutation(n, g);

  sigma_w.change(i - 1, sigma_w.inverse()[sigma[a]]);
  sigma_w.change(i, sigma_w.inverse()[sigma[b]]);

  return sigma_w;
}
} // namespace internal

//! Generates cloaking element of the form w b_i^{\pm 2} w^{-1} for a permutation p.
//! The word w is of random length [min_length, max_length].
//! Numbers 0 \leq a < b \leq n-1 are relates to t-values used from evaluation of E-multiplication.
template <typename URNG>
Word generateCloakingElement(
    size_t n,
    size_t a,
    size_t b,
    size_t min_length,
    size_t max_length,
    const Permutation& p,
    int center_degree,
    URNG& g) {
  assert(n > 3);

  assert(0 <= a);
  assert(a < b);
  assert((b + 1) <= n);

  assert(min_length <= max_length);

  const auto abs_center_degree = std::abs(center_degree);

  assert(abs_center_degree == 2 || abs_center_degree == 4);

  boost::random::uniform_int_distribution<size_t> dist(1, n - 1);

  return internal::generateCloakingElement(n, p[a], p[b], dist(g), min_length, max_length, center_degree, g);
}

//! Generates a product of cloaking_number cloaking elements.
template <typename URNG>
Word generateCloakingElement(
    size_t cloaking_number,
    size_t n,
    size_t a,
    size_t b,
    size_t min_length,
    size_t max_length,
    const Permutation& p,
    int center_degree,
    URNG& g) {
  Word result;

  for (size_t i = 0; i < cloaking_number; ++i) {
    result *= generateCloakingElement(n, a, b, min_length, max_length, p, center_degree, g);
  }

  return result;
}

//! Canonical cloaking element for sigma w.r.t. Walnut spec.
template <typename URNG>
Word generateCanonicalCloakingEl(
    size_t n, size_t a, size_t b, size_t L, const Permutation& sigma, size_t center_degree, URNG& g) {
  assert(n > 3);

  assert(0 <= a);
  assert(a < b);
  assert((b + 1) <= n);

  assert((center_degree == 2) || (center_degree == 4));

  boost::random::uniform_int_distribution<int> dist(2, n - 1);
  const auto i = dist(g);

  const auto sigma_w = internal::generateSigmaW(n, a, b, i, sigma, g);
  auto w = -random::randomWord(sigma_w, g);

  // append L random pure braid subgroup generators
  for (size_t j = 0; j < L; ++j) {
    const auto ij = random::subset<size_t>(2, 1, n, g);
    auto pure_braid = braidgroup::getPureBraidSubgroupGenerator(n, *ij.begin(), *(++ij.begin()));

    if (g() % 2 == 0) {
      pure_braid = -pure_braid;
    }

    w *= pure_braid;
  }

  Word center = Word(i) ^ center_degree;

  if (g() % 2 != 0) {
    center = -center;
  }

  return w * center * -w;
}
} // namespace coloredburau
} // namespace crag

#endif // CRAG_CLOAKING_ELEMENT_H
