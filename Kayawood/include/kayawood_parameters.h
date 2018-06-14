#pragma once

#ifndef CRAG_KAYAWOOD_PARAMETERS_H
#define CRAG_KAYAWOOD_PARAMETERS_H

#include <boost/random/uniform_int_distribution.hpp>
#include <vector>

#include "FiniteField.h"
#include "random_subset.h"

namespace crag {
namespace kayawood {

//! Represents public parameters of the Kayawood protocol:
//! 1. n - rank of braid group B_n.
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
      , r_(10)
      , z_min_length_(300)
      , z_max_length_(400)
      , beta_min_length_(10)
      , beta_max_length_(15)
      , alpha_min_length_(300)
      , alpha_max_length_(400)
      , bob_private_key_min_length_(22)
      , bob_private_key_max_length_(22) {
    if ((n_ < 16) || (n_ % 2 == 1)) {
      throw std::invalid_argument("Expect n to be an even number >= 16.");
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

  size_t r() const {
    return r_;
  }

  size_t zMinLength() const {
    return z_min_length_;
  }

  size_t zMaxLength() const {
    return z_max_length_;
  }

  size_t betaMinLength() const {
    return beta_min_length_;
  }

  size_t betaMaxLength() const {
    return beta_max_length_;
  }

  size_t alphaMinLength() const {
    return alpha_min_length_;
  }

  size_t alphaMaxLength() const {
    return alpha_max_length_;
  }

  size_t bobPrivateKeyMinLength() const {
    return bob_private_key_min_length_;
  }

  size_t bobPrivateKeyMaxLength() const {
    return bob_private_key_max_length_;
  }

  PublicParameters& r(size_t r) {
    r_ = r;
    return *this;
  }

  PublicParameters& zMinLength(size_t min_length) {
    z_min_length_ = min_length;
    z_max_length_ = std::max(z_max_length_, z_min_length_);

    return *this;
  }

  PublicParameters& zMaxLength(size_t max_length) {
    z_max_length_ = max_length;
    z_min_length_ = std::min(z_min_length_, z_max_length_);

    return *this;
  }

  PublicParameters& betaMinLength(size_t min_length) {
    beta_min_length_ = min_length;
    beta_max_length_ = std::max(beta_max_length_, beta_min_length_);

    return *this;
  }

  PublicParameters& betaMaxLength(size_t max_length) {
    beta_max_length_ = max_length;
    beta_min_length_ = std::min(beta_min_length_, beta_max_length_);

    return *this;
  }

  PublicParameters& alphaMinLength(size_t min_length) {
    alpha_min_length_ = min_length;
    alpha_max_length_ = std::max(alpha_max_length_, alpha_min_length_);

    return *this;
  }

  PublicParameters& alphaMaxLength(size_t max_length) {
    alpha_max_length_ = max_length;
    alpha_min_length_ = std::min(alpha_min_length_, alpha_max_length_);

    return *this;
  }

  PublicParameters& bobPrivateKeyMinLength(size_t min_length) {
    bob_private_key_min_length_ = min_length;
    bob_private_key_max_length_ = std::max(bob_private_key_max_length_, bob_private_key_min_length_);

    return *this;
  }

  PublicParameters& bobPrivateKeyMaxLength(size_t max_length) {
    bob_private_key_max_length_ = max_length;
    bob_private_key_min_length_ = std::min(bob_private_key_min_length_, bob_private_key_max_length_);

    return *this;
  }

private:
  size_t n_;
  size_t a_;
  size_t b_;

  std::vector<T> t_values_;

  size_t r_;

  size_t z_min_length_;
  size_t z_max_length_;

  size_t beta_min_length_;
  size_t beta_max_length_;

  size_t alpha_min_length_;
  size_t alpha_max_length_;

  size_t bob_private_key_min_length_;
  size_t bob_private_key_max_length_;
};

template <typename T, typename Stabilizer>
class TValuesGenerator {};

class StabilizerSquare;
class StabilizerManyShort;
class StabilizerDoubleSquare;
class StabilizerTrivial;

template <typename T>
class TValuesGenerator<T, StabilizerSquare> {
public:
  TValuesGenerator(size_t n, size_t a, size_t b)
      : n_(n)
      , a_(a)
      , b_(b) {}

  template <typename URNG>
  std::vector<T> operator()(URNG& g) const {
    std::vector<T> t_values(n_, T(0));

    t_values[a_] = t_values[b_] = T(1);

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
class TValuesGenerator<T, StabilizerManyShort> {
public:
  TValuesGenerator(size_t n, size_t a, size_t b)
      : n_(n)
      , a_(a)
      , b_(b) {}

  template <typename URNG>
  std::vector<T> operator()(URNG& g) const {
    return TValuesGenerator<T, StabilizerSquare>(n_, a_, b_)(g);
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
    std::vector<T> t_values(n_, T(0));

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

template <typename T>
class TValuesGenerator<T, StabilizerTrivial> {
public:
  TValuesGenerator(size_t n, size_t a, size_t b)
      : n_(n) {}

  template <typename URNG>
  std::vector<T> operator()(URNG& g) const {
    std::vector<T> t_values(n_, T(0));

    for (size_t i = 0; i < n_; ++i) {
      t_values[i] = finitefield::generateNonZeroNonUnit<T>(g);
    }

    return t_values;
  }

private:
  size_t n_;
};

//! Given the rank n of the Braid group, generates random public parameters for Kayawood protocol.
template <typename T, typename Stabilizer, typename URNG>
PublicParameters<T> randomParameters(size_t n, URNG& g) {
  const auto a_b_set = random::subset<size_t>(2, 0, n - 1, g);
  const auto a = *a_b_set.begin();
  const auto b = *a_b_set.rbegin();

  auto t_values = TValuesGenerator<T, Stabilizer>(n, a, b)(g);

  return PublicParameters<T>(n, a, b, std::move(t_values));
}

template <typename T, typename Stabilizer>
PublicParameters<T> randomParameters(size_t n, size_t seed) {
  std::mt19937_64 g(seed);
  return randomParameters<T, Stabilizer>(n, g);
}
} // namespace kayawood
} // namespace crag

#endif // CRAG_KAYAWOOD_PARAMETERS_H
