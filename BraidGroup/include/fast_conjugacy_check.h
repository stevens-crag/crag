#pragma once

#ifndef CRAG_FAST_CONJUGACY_CHECK_H
#define CRAG_FAST_CONJUGACY_CHECK_H

#include <boost/functional/hash.hpp>

#include "colored_burau.h"
#include "matrix.h"

namespace crag {
namespace braidgroup {

template <typename T>
class FastConjugacyChecker {
public:
  explicit FastConjugacyChecker(size_t n)
      : unit_(getUnitEl_(n, 0)) {}

  FastConjugacyChecker(size_t n, size_t seed)
      : unit_(getUnitEl_(n, seed)) {}

  template <typename URNG>
  FastConjugacyChecker(size_t n, URNG& g)
      : unit_(getUnitEl_(n, g)) {}

  //! Returns true iff lhs and rhs are not conjugate.
  //! Returning false means that lhs and rhs are most likely conjugate.
  bool areNotConjugate(const Word& lhs, const Word& rhs) const {
    return matrix::charPoly((unit_ * lhs).matrix()) != matrix::charPoly((unit_ * rhs).matrix());
  }

private:
  coloredburau::CBProjectionElement<T> unit_;

  template <typename URNG>
  coloredburau::CBProjectionElement<T> getUnitEl_(size_t n, URNG& g) const {
    if (n < 3) {
      throw std::invalid_argument("Expect n to be greater or equal to 3.");
    }

    std::vector<T> t_values(n, finitefield::generateNonZeroNonUnit<T>(g));

    return coloredburau::CBProjectionElement<T>(std::move(t_values));
  }

  coloredburau::CBProjectionElement<T> getUnitEl_(size_t n, size_t seed) const {
    std::mt19937_64 g(seed);
    return getUnitEl_(n, g);
  }
};
} // namespace braidgroup
} // namespace crag

#endif // CRAG_FAST_CONJUGACY_CHECK_H
