#pragma once

#ifndef CRAG_MONOMIAL_H
#define CRAG_MONOMIAL_H

#include <algorithm>
#include <type_traits>
#include <vector>

#include "FiniteField.h"

namespace crag {
namespace polynomials {

template <typename ExponentType>
using Term = std::vector<ExponentType>;

template <typename T, typename ExponentType>
class Monomial {
public:
  using exponent_t = ExponentType;
  using term_t = Term<ExponentType>;

  Monomial() = delete;

  explicit Monomial(term_t term)
      : coef_(1)
      , term_(std::move(term)) {
    validateTerm_();
  }

  Monomial(T coef, term_t term)
      : coef_(std::move(coef))
      , term_(std::move(term)) {
    validateTerm_();
  }

  const T& coef() const {
    return coef_;
  }

  const term_t& term() const {
    return term_;
  }

  size_t dimension() const {
    return term_.size();
  }

  bool isZero() const {
    return coef_ == 0;
  }

  //! Returns true if all variables' exponents are zeros
  bool isNumber() const {
    return isZero() || std::all_of(term_.begin(), term_.end(), [](const ExponentType& e) { return e == 0; });
  }

  bool isUnit() const {
    return (coef_ == 1) && isNumber();
  }

  Monomial& operator*=(const term_t& term) {
    if (dimension() != term.size()) {
      throw std::invalid_argument("Dimensions of terms don't match.");
    }

    for (size_t i = 0; i < term_.size(); ++i) {
      term_[i] += term[i];
    }

    return *this;
  }

  Monomial& operator*=(const T& coef) {
    coef_ *= coef;
    return *this;
  }

  Monomial& operator*=(const Monomial& other) {
    return ((*this) *= other.term()) *= other.coef();
  }

  Monomial operator-() const {
    return Monomial(-coef_, term_);
  }

  //! Laurent monomials can be divided by terms
  template <
      typename S, typename = typename std::enable_if<std::is_same<S, ExponentType>::value>::type,
      typename = typename std::enable_if<std::is_signed<S>::value>::type>
  Monomial& operator/=(const Term<S>& term) {
    auto inverse_term = term;

    for (auto& exponent : inverse_term) {
      exponent *= -1;
    }

    return (*this) *= inverse_term;
  }

  bool operator==(const Monomial& other) const {
    if (dimension() != other.dimension()) {
      throw std::invalid_argument("Dimensions of monomials don't match.");
    }

    if (other.isZero()) {
      return isZero();
    }

    return (coef_ == other.coef()) && (term_ == other.term());
  }

  bool operator!=(const Monomial& other) const {
    return !(*this == other);
  }

  bool operator==(const term_t& term) const {
    if (dimension() != term.size()) {
      throw std::invalid_argument("Dimensions of monomials don't match.");
    }

    return (coef_ == 1) && (term_ == term);
  }

  bool operator!=(const term_t& term) const {
    return !(*this == term);
  }

  bool operator==(const T& coef) const {
    if (coef == 0) {
      return isZero();
    }

    return (coef_ == coef) && isNumber();
  }

  bool operator!=(const T& coef) const {
    return !(*this == coef);
  }

private:
  T coef_;
  term_t term_;

  void validateTerm_() const {
    if (term_.size() < 2) {
      throw std::invalid_argument("Terms must have at least 2 variables.");
    }
  }
};

template <
    typename T, typename ExponentType, typename = typename std::enable_if<std::is_signed<ExponentType>::value>::type>
Monomial<T, ExponentType> operator/(Monomial<T, ExponentType> monomial, const Term<ExponentType>& term) {
  monomial /= term;
  return monomial;
}

template <typename T, typename ExponentType>
Monomial<T, ExponentType> operator*(Monomial<T, ExponentType> lhs, const Monomial<T, ExponentType>& rhs) {
  lhs *= rhs;
  return lhs;
}

template <typename T, typename ExponentType>
Monomial<T, ExponentType> operator*(Monomial<T, ExponentType> lhs, const Term<ExponentType>& term) {
  lhs *= term;
  return lhs;
}

template <typename T, typename ExponentType>
Monomial<T, ExponentType> operator*(const Term<ExponentType>& term, Monomial<T, ExponentType> monomial) {
  monomial *= term;
  return monomial;
}

template <typename T, typename ExponentType>
Monomial<T, ExponentType> operator*(Monomial<T, ExponentType> lhs, const T& coef) {
  lhs *= coef;
  return lhs;
}

template <typename T, typename ExponentType>
Monomial<T, ExponentType> operator*(const T& coef, Monomial<T, ExponentType> monomial) {
  monomial *= coef;
  return monomial;
}

template <typename T, typename ExponentType>
bool operator==(const Term<ExponentType>& term, const Monomial<T, ExponentType>& monomial) {
  return monomial == term;
}

template <typename T, typename ExponentType>
bool operator!=(const Term<ExponentType>& term, const Monomial<T, ExponentType>& monomial) {
  return !(monomial == term);
}

template <typename T, typename ExponentType>
bool operator==(const T& coef, const Monomial<T, ExponentType>& monomial) {
  return monomial == coef;
}

template <typename T, typename ExponentType>
bool operator!=(const T& coef, const Monomial<T, ExponentType>& monomial) {
  return !(monomial == coef);
}

//! Evaluates a monomial at a point.
//! May throw if monomial's term has negative exponents and the corresponding value is zero.
template <typename T, typename ExponentType>
T evaluate(const Monomial<T, ExponentType>& monomial, const std::vector<T>& values);

template <typename T, typename ExponentType>
T evaluate(const T& coef, const Term<ExponentType>& term, const std::vector<T>& values);

using finitefield::FieldElement;

//! Implementation for finite fields
template <typename Ideal, typename ExponentType>
FieldElement<Ideal> evaluate(
    const FieldElement<Ideal>& coef, const Term<ExponentType>& term, const std::vector<FieldElement<Ideal>>& values) {
  if (term.size() != values.size()) {
    throw std::invalid_argument("The number of values provided doesn't match dimension of the term.");
  }

  if (coef == FieldElement<Ideal>(0)) {
    return coef;
  }

  FieldElement<Ideal> result(coef);

  for (size_t i = 0; i < term.size(); ++i) {
    result *= finitefield::pwr(values[i], term[i]);
  }

  return result;
}

template <typename Ideal, typename ExponentType>
FieldElement<Ideal>
evaluate(const Monomial<FieldElement<Ideal>, ExponentType>& monomial, const std::vector<FieldElement<Ideal>>& values) {
  return evaluate(monomial.coef(), monomial.term(), values);
}

//! Monomial of standard polynomial
template <typename T>
using StandardMonomial = Monomial<T, uint32_t>;
using IntStandardMonomial = StandardMonomial<int>;

//! Laurent monomial
template <typename T>
using LaurentMonomial = Monomial<T, int32_t>;
using IntLaurentMonomial = LaurentMonomial<int>;
} // namespace polynomials
} // namespace crag

#endif // CRAG_MONOMIAL_H
