#pragma once

#ifndef CRAG_WALNUT_CLOAKING_H
#define CRAG_WALNUT_CLOAKING_H

#include "Permutation.h"
#include "Word.h"
#include "cloaking_element.h"
#include "walnut_parameters.h"

namespace crag {
namespace walnut {

//! Generator of elements from the stabilizer w.r.t E-multiplication.
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
  Word operator()(const Permutation& p, URNG& g) const {
    return coloredburau::generateCanonicalCloakingEl(n_, a_, b_, L_, p, 2, g);
  }

private:
  size_t n_;
  size_t a_;
  size_t b_;
  size_t L_;
};

//! Generator of elements from the stabilizer w.r.t E-multiplication.
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
  Word operator()(const Permutation& p, URNG& g) const {
    return coloredburau::generateCanonicalCloakingEl(n_, a_, b_, L_, p, 4, g);
  }

private:
  size_t n_;
  size_t a_;
  size_t b_;
  size_t L_;
};
} // namespace walnut
} // namespace crag

#endif // CRAG_WALNUT_CLOAKING_H
