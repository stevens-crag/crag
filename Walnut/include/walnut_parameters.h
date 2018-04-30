#pragma once

#ifndef CRAG_WALNUT_PARAMETERS_H
#define CRAG_WALNUT_PARAMETERS_H

#include <vector>

#include "FiniteField.h"
#include "random_subset.h"

namespace crag {
namespace walnut {

//! Represents public parameters of the Walnut protocol:
//! 1. n - rank of braid group B_n, n >= 8.
//! 2. Two numbers 0 \leq a < b \leq n - 1.
//! 3. T-values {\tau_0,\dots,\tau_{n-1}} \in T^n,
//!     properties of \tau_a, \tau_b depends on algorithms used for generation of cloaking elements.
template <typename T>
class PublicParameters {
public:
  PublicParameters(size_t n, size_t a, size_t b, std::vector<T> t_values)
      : n_(n)
      , a_(a)
      , b_(b)
      , t_values_(std::move(t_values))
      , w_min_length_(280)
      , w_max_length_(300) {
    if (n_ < 8) {
      throw std::invalid_argument("Expect n to be an even number >= 8.");
    }

    if (!((0 <= a_) && (a_ < b_) && (b_ + 1 <= n_))) {
      throw std::invalid_argument("Invalid values a and/or b.");
    }

    if (t_values_.size() != n_) {
      throw std::invalid_argument("Number of t-values must be equal to n.");
    }
  }

  size_t n() const {
    return n_;
  }

  size_t a() const {
    return a_;
  }

  size_t b() const {
    return b_;
  }

  const std::vector<T>& tValues() const {
    return t_values_;
  }

  size_t wMinLength() const {
    return w_min_length_;
  }

  size_t wMaxLength() const {
    return w_max_length_;
  }

  PublicParameters& wMinLength(size_t min_length) {
    w_min_length_ = min_length;
    w_max_length_ = std::max(w_max_length_, w_min_length_);

    return *this;
  }

  PublicParameters& wMaxLength(size_t max_length) {
    w_max_length_ = max_length;
    w_min_length_ = std::min(w_min_length_, w_max_length_);

    return *this;
  }

private:
  size_t n_;
  size_t a_;
  size_t b_;

  std::vector<T> t_values_;

  size_t w_min_length_;
  size_t w_max_length_;
};

template <typename T, typename Stabilizer>
class TValuesGenerator {};

class StabilizerSquare;
class StabilizerDoubleSquare;

template <typename T>
class TValuesGenerator<T, StabilizerSquare> {
public:
  TValuesGenerator(size_t n, size_t a, size_t b)
      : n_(n)
      , a_(a)
      , b_(b) {}

  template <typename URNG>
  std::vector<T> operator()(URNG& g) const {
    const T zero(0);
    const T unit(1);

    std::vector<T> t_values(n_, zero);

    t_values[a_] = t_values[b_] = unit;

    for (size_t i = 0; i < n_; ++i) {
      if ((i == a_) || (i == b_)) {
        continue;
      }

      t_values[i] = finitefield::generateNonZeroNonUnit<T>(g);
    }

    return t_values;
  }

private:
  size_t n_;
  size_t a_;
  size_t b_;
};

template <typename T>
class TValuesGenerator<T, StabilizerDoubleSquare> {
public:
  TValuesGenerator(size_t n, size_t a, size_t b)
      : n_(n)
      , a_(a)
      , b_(b) {}

  template <typename URNG>
  std::vector<T> operator()(URNG& g) const {
    const T zero(0);
    const T unit(1);

    std::vector<T> t_values(n_, zero);

    t_values[a_] = finitefield::generateNonZeroNonUnit<T>(g);
    t_values[b_] = -t_values[a_].inverse();

    for (size_t i = 0; i < n_; ++i) {
      if ((i == a_) || (i == b_)) {
        continue;
      }

      t_values[i] = finitefield::generateNonZeroNonUnit<T>(g);
    }

    return t_values;
  }

private:
  size_t n_;
  size_t a_;
  size_t b_;
};

//! Given the rank n of the Braid group, generates random public parameters for Walnut protocol.
template <typename T, typename Stabilizer, typename URNG>
PublicParameters<T> randomParameters(size_t n, URNG& g) {
  if (n < 8) {
    throw std::invalid_argument("Expect n to be an even number >= 8.");
  }

  const auto a_b_set = random::subset<size_t>(2, 0, n - 1, g);
  const auto a = *a_b_set.begin();
  const auto b = *a_b_set.rbegin();

  const T zero(0);
  const T unit(1);

  auto t_values = TValuesGenerator<T, Stabilizer>(n, a, b)(g);

  return PublicParameters<T>(n, a, b, std::move(t_values));
}

template <typename T, typename Stabilizer>
PublicParameters<T> randomParameters(size_t n, size_t seed) {
  std::mt19937_64 g(seed);
  return randomParameters<T, Stabilizer>(n, g);
}
} // namespace walnut
} // namespace crag

#endif // CRAG_WALNUT_PARAMETERS_H
