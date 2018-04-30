#pragma once

#ifndef CRAG_FAST_IDENTITY_CHECK_H
#define CRAG_FAST_IDENTITY_CHECK_H

#include <boost/functional/hash.hpp>

#include "colored_burau.h"

namespace crag {
namespace braidgroup {

//! Fast identity checker that uses colored Burau mapping over T.
template <typename T>
class FastIdentityChecker {
public:
  explicit FastIdentityChecker(size_t n)
      : unit_(getUnitEl_(n, 0)) {}

  FastIdentityChecker(size_t n, size_t seed)
      : unit_(getUnitEl_(n, seed)) {}

  template <typename URNG>
  FastIdentityChecker(size_t n, URNG& g)
      : unit_(getUnitEl_(n, g)) {}

  //! Applies E-multiplication by w to the trivial pair (E, id)
  //! and returns true iff the result is not the trivial pair,
  //! which means that w is not a trivial braid.
  bool isNonTrivial(const Word& w) const {
    return unit_ != (unit_ * w);
  }

private:
  coloredburau::CBProjectionElement<T> unit_;

  template <typename URNG>
  coloredburau::CBProjectionElement<T> getUnitEl_(size_t n, URNG& g) const {
    if (n < 3) {
      throw std::invalid_argument("Expect n to be greater or equal to 3.");
    }

    std::vector<T> t_values;
    t_values.reserve(n);

    for (size_t i = 0; i < n; ++i) {
      t_values.push_back(finitefield::generateNonZeroNonUnit<T>(g));
    }

    return coloredburau::CBProjectionElement<T>(std::move(t_values));
  }

  coloredburau::CBProjectionElement<T> getUnitEl_(size_t n, size_t seed) const {
    std::mt19937_64 g(seed);
    return getUnitEl_(n, g);
  }
};

//! Computes the hash of CBProjectionElement<T> corresponding to a word w.
template <typename T>
size_t projectionHash(const Word& w, const std::vector<T>& t_values) {
  return std::hash<coloredburau::CBProjectionElement<T>>()(coloredburau::project(w, t_values));
}

template <typename T>
class BraidHasher {
public:
  using braid_hash_t = size_t;

  explicit BraidHasher(size_t n)
      : n_(n)
      , t_values(generateTValues_(n, 0)) {}

  BraidHasher(size_t n, size_t seed)
      : n_(n)
      , t_values(generateTValues_(n, seed)) {}

  template <typename URNG>
  BraidHasher(size_t n, URNG& g)
      : n_(n)
      , t_values(generateTValues_(n, g)) {}

  braid_hash_t operator()(const Word& w) const {
    return projectionHash(w, t_values);
  }

  braid_hash_t operator()(const std::vector<Word>& words) const {
    std::vector<braid_hash_t> hashes;
    hashes.reserve(words.size());

    for (const auto& w : words) {
      hashes.push_back((*this)(w));
    }

    return boost::hash_value(hashes);
  }

private:
  size_t n_;
  std::vector<T> t_values;

  template <typename URNG>
  std::vector<T> generateTValues_(size_t n, URNG& g) const {
    std::vector<T> t_values;
    t_values.reserve(n);

    for (size_t i = 0; i < n; ++i) {
      t_values.push_back(finitefield::generateNonZeroNonUnit<T>(g));
    }

    return t_values;
  }

  std::vector<T> generateTValues_(size_t n, size_t seed) const {
    std::mt19937_64 g(seed);
    return generateTValues_(n, g);
  }
};
} // namespace braidgroup
} // namespace crag

#endif // CRAG_FAST_IDENTITY_CHECK_H
