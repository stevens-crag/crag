#pragma once

#ifndef CRAG_GCD_H
#define CRAG_GCD_H

#include <boost/math/tools/polynomial.hpp>

namespace crag {
namespace finitefield {

//! Returns pair a and b, where a * b = x, b is invertible and a is "normed".
//! For example, if x in ZZ, then a is the abs of x, and b is the sign of x.
//! If x in F[x], then a is the corresponding unitary polynomial, and b is the lead coefficient.
template <typename RingElement>
std::pair<RingElement, RingElement> getNormed(const RingElement& x);

template <>
std::pair<int, int> getNormed<int>(const int& x);

template <typename T>
std::pair<boost::math::tools::polynomial<T>, boost::math::tools::polynomial<T>>
getNormed(const boost::math::tools::polynomial<T>& x) {
  const auto lead = x[x.degree()];
  return std::make_pair(x * (T(1) / lead), boost::math::tools::polynomial<T>(lead));
}


//! Returns gcd d and x, y s. t. d = ax + by.
template <typename RingElement>
RingElement getGCD(const RingElement& a, const RingElement& b, RingElement& x, RingElement& y) {
  if (b == RingElement(0)) {
    x = RingElement(1);
    y = RingElement(0);

    return a;
  }

  const auto d = getGCD(b, a % b, x, y);
  const auto t = x;
  x = y;
  y = t - (a / b) * y;

  return d;
}

} // namespace finitefield
} // namespace crag

#endif /* CRAG_GCD_H */
