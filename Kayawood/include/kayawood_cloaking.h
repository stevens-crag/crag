#pragma once

#ifndef CRAG_KAYAWOOD_CLOAKING_H
#define CRAG_KAYAWOOD_CLOAKING_H

#include <boost/random/uniform_int_distribution.hpp>

#include "Permutation.h"
#include "Word.h"
#include "cloaking_element.h"
#include "kayawood_parameters.h"

namespace crag {
namespace kayawood {

// Stabilizer objects take \sigma and w and return u such that (M, \sigma) * w = (M, \sigma) * u

//! Generator of elements from the stabilizer w.r.t. E-multiplication.
//! Generates elements of the form w x_i^{+/- 2} w^{-1}
class StabilizerSquare {
public:
  template <typename T>
  StabilizerSquare(const PublicParameters<T>& parameters, size_t L)
      : n_(parameters.n())
      , a_(parameters.a())
      , b_(parameters.b())
      , L_(L) {
    if ((parameters.tValues()[a_] != T(1)) || (parameters.tValues()[b_] != T(1))) {
      throw std::invalid_argument("T-values in positions a and b must be equal to 1.");
    }
  }

  template <typename URNG>
  Word operator()(const Permutation& p, const Word& w, URNG& g) const {
    const auto left_el = coloredburau::generateCanonicalCloakingEl(n_, a_, b_, L_, p, 2, g);

    const auto right_p = coloredburau::permutation(n_, w) * p;
    const auto right_el = coloredburau::generateCanonicalCloakingEl(n_, a_, b_, L_, right_p, 2, g);

    return left_el * w * right_el;
  }

private:
  size_t n_;
  size_t a_;
  size_t b_;
  size_t L_;
};

//! Generator of elements from the stabilizer w.r.t. E-multiplication.
//! Generates elements of the form \prod w x_i^{+/- 2} w^{-1}.
//! Length of w and number of factors can be specified.
class StabilizerManyShort {
public:
  template <typename T>
  StabilizerManyShort(const PublicParameters<T>& parameters, size_t w_min_length, size_t w_max_length, size_t count)
      : n_(parameters.n())
      , a_(parameters.a())
      , b_(parameters.b())
      , w_min_length_(w_min_length)
      , w_max_length_(w_max_length)
      , count_(count) {
    if ((parameters.tValues()[a_] != T(1)) || (parameters.tValues()[b_] != T(1))) {
      throw std::invalid_argument("T-values in positions a and b must be equal to 1.");
    }
  }

  template <typename URNG>
  Word operator()(const Permutation& p, const Word& w, URNG& g) const {
    Word result = w;

    for (size_t i = 0; i < count_; ++i) {
      result = insertCloakingElAtRandomPosition_(p, result, g);
    }

    return result;
  }

private:
  size_t n_;
  size_t a_;
  size_t b_;
  size_t w_min_length_;
  size_t w_max_length_;
  size_t count_;

  template <typename URNG>
  Word insertCloakingElAtRandomPosition_(const Permutation& p, const Word& w, URNG& g) const {
    assert(!w.empty());

    const auto distribution = boost::random::uniform_int_distribution<size_t>(0, w.length());
    const auto position = distribution(g);

    const auto w_left = w.subword(0, position);
    const auto w_right = w.subword(position, w.length());

    const auto sigma = p * coloredburau::permutation(n_, w_left);

    const auto cloaking_el =
        coloredburau::generateCloakingElement(1, n_, a_, b_, w_min_length_, w_max_length_, sigma, 2, g);

    return w_left * cloaking_el * w_right;
  }
};

//! Generator of elements from the stabilizer w.r.t. E-multiplication.
//! Generates elements of the form w x_i^{+/- 4} w^{-1}
class StabilizerDoubleSquare {
public:
  template <typename T>
  StabilizerDoubleSquare(const PublicParameters<T>& parameters, size_t L)
      : n_(parameters.n())
      , a_(parameters.a())
      , b_(parameters.b())
      , L_(L) {
    if (parameters.tValues()[a_] * parameters.tValues()[b_] != -T(1)) {
      throw std::invalid_argument("T-values in positions a and b must give -1 in product.");
    }
  }

  template <typename URNG>
  Word operator()(const Permutation& p, const Word& w, URNG& g) const {
    const auto left_el = coloredburau::generateCanonicalCloakingEl(n_, a_, b_, L_, p, 4, g);

    const auto right_p = coloredburau::permutation(n_, w) * p;
    const auto right_el = coloredburau::generateCanonicalCloakingEl(n_, a_, b_, L_, right_p, 4, g);

    return left_el * w * right_el;
  }

private:
  size_t n_;
  size_t a_;
  size_t b_;
  size_t L_;
};

class StabilizerTrivial {
public:
  template <typename URNG>
  Word operator()(const Permutation& p, const Word& w, URNG& g) const {
    return w;
  }
};
} // namespace kayawood
} // namespace crag

#endif // CRAG_KAYAWOOD_CLOAKING_H
