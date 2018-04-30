#pragma once

#ifndef CRAG_POLYNOMIAL_H
#define CRAG_POLYNOMIAL_H

#include <functional>
#include <map>
#include <sstream>
#include <string>
#include <type_traits>
#include <unordered_set>
#include <vector>

#include "monomial.h"

namespace crag {
namespace polynomials {

template <typename T, typename ExponentType, typename TermCompare = std::less<Term<ExponentType>>>
class MonomialSum {
public:
  using exponent_t = ExponentType;
  using term_t = Term<ExponentType>;
  using monomial_t = Monomial<T, ExponentType>;
  using MonomialMap = std::map<term_t, T, TermCompare>;
  using const_iterator = typename MonomialMap::const_iterator;

  //! Constructs empty monomial sum of given dimension
  explicit MonomialSum(size_t dimension)
      : dimension_(dimension) {
    if (dimension_ < 2) {
      throw std::invalid_argument("Multivariate polynomial requires at least 2 variables.");
    }
  }

  explicit MonomialSum(const monomial_t& monomial)
      : dimension_(monomial.dimension())
      , monomials_({{monomial.term(), monomial.coef()}}) {}

  MonomialSum(const T& coef, const term_t& term)
      : dimension_(term.size())
      , monomials_({{term, coef}}) {
    if (dimension_ < 2) {
      throw std::invalid_argument("Multivariate polynomial requires at least 2 variables.");
    }
  }

  //! Returns the number of monomials
  size_t size() const {
    return monomials_.size();
  }

  //! Returns the number of variables
  size_t dimension() const {
    return dimension_;
  }

  const_iterator begin() const {
    return monomials_.begin();
  }

  const_iterator end() const {
    return monomials_.end();
  }

  //! Checks if the sum is actually zero
  bool isZero() const {
    return monomials_.empty();
  }

  //! Checks if the sum is actually a number
  bool isNumber() const {
    if (isZero()) {
      return true;
    }

    if (size() != 1) {
      return false;
    }

    const auto term = monomials_.begin()->first;

    return std::all_of(term.begin(), term.end(), [](const ExponentType& e) { return e == 0; });
  }

  //! Checks if the sum is actually unit
  bool isUnit() const {
    return (monomials_.size() == 1) && (monomials_.begin()->second == T(1)) && isNumber();
  }

  MonomialSum& operator+=(const monomial_t& monomial) {
    if (dimension() != monomial.dimension()) {
      throw std::invalid_argument("Dimensions of polynomials don't match.");
    }

    if (monomial.coef() == T(0)) {
      return *this;
    }

    auto it = monomials_.find(monomial.term());

    if (it == monomials_.end()) {
      // add new term
      monomials_.emplace(monomial.term(), monomial.coef());
    } else {
      // update coef of existing term
      it->second += monomial.coef();

      if (it->second == T(0)) {
        monomials_.erase(it);
      }
    }

    return *this;
  }

  MonomialSum& operator+=(const MonomialSum& other) {
    if (dimension_ != other.dimension_) {
      throw std::invalid_argument("Dimensions of polynomials don't match.");
    }

    for (const auto& it : other) {
      (*this) += monomial_t(it.second, it.first);
    }

    return *this;
  }

  MonomialSum& operator+=(const term_t& term) {
    return (*this) += monomial_t(T(1), term);
  }

  MonomialSum& operator+=(const T& coef) {
    return (*this) += monomial_t(coef, term_t(dimension(), 0));
  }

  MonomialSum& operator-=(const monomial_t& monomial) {
    return (*this) += -monomial;
  }

  MonomialSum& operator-=(const MonomialSum& other) {
    if (dimension_ != other.dimension_) {
      throw std::invalid_argument("Dimensions of polynomials don't match.");
    }

    for (const auto& it : other) {
      (*this) += monomial_t(-it.second, it.first);
    }

    return *this;
  }

  MonomialSum& operator-=(const term_t& term) {
    return (*this) += monomial_t(-T(1), term);
  }

  MonomialSum& operator-=(const T& coef) {
    return (*this) += monomial_t(-coef, term_t(dimension(), 0));
  }

  MonomialSum operator-() const {
    auto result = *this;
    return result *= -T(1);
  }

  MonomialSum& operator*=(const MonomialSum& other) {
    if (dimension_ != other.dimension_) {
      throw std::invalid_argument("Dimensions of polynomials don't match.");
    }

    MonomialSum result(dimension_);

    for (const auto& this_monomial : monomials_) {
      for (const auto& other_monomial : other.monomials_) {
        auto m = monomial_t(this_monomial.second, this_monomial.first);

        m *= other_monomial.first;
        m *= other_monomial.second;

        result += m;
      }
    }

    monomials_ = result.monomials_;

    return *this;
  }

  MonomialSum& operator*=(const monomial_t& monomial) {
    if (dimension_ != monomial.dimension()) {
      throw std::invalid_argument("Dimensions of polynomials don't match.");
    }

    return (*this) *= (MonomialSum(dimension_) + monomial);
  }

  MonomialSum& operator*=(const term_t& term) {
    return (*this) *= monomial_t(T(1), term);
  }

  MonomialSum& operator*=(const T& coef) {
    if (coef == T(0)) {
      monomials_.clear();
    }

    for (auto& it : monomials_) {
      it.second *= coef;
    }

    return *this;
  }

  //! Laurent polynomials can be divided by terms
  template <
      typename S, typename = typename std::enable_if<std::is_same<S, ExponentType>::value>::type,
      typename = typename std::enable_if<std::is_signed<S>::value>::type>
  MonomialSum& operator/=(const Term<S>& term) {
    auto inverse_term = term;

    for (auto& exponent : inverse_term) {
      exponent *= -1;
    }

    return (*this) *= inverse_term;
  }

  bool operator==(const MonomialSum& other) const {
    if (dimension_ != other.dimension_) {
      throw std::invalid_argument("Dimensions of polynomials don't match.");
    }

    return monomials_ == other.monomials_;
  }

  bool operator!=(const MonomialSum& other) const {
    return !(*this == other);
  }

  bool operator==(const monomial_t& monomial) const {
    if (dimension() != monomial.dimension()) {
      throw std::invalid_argument("Dimensions of polynomials don't match.");
    }

    if (monomial.isZero()) {
      return isZero();
    }

    if (size() != 1) {
      return false;
    }

    const auto it = monomials_.begin();

    return (it->second == monomial.coef()) && (it->first == monomial.term());
  }

  bool operator!=(const monomial_t& monomial) const {
    return !(*this == monomial);
  }

  bool operator==(const term_t& term) const {
    return *this == monomial_t(T(1), term);
  }

  bool operator!=(const term_t& term) const {
    return !(*this == monomial_t(T(1), term));
  }

  bool operator==(const T& coef) const {
    if (coef == T(0)) {
      return isZero();
    }

    return isNumber() && (monomials_.begin()->second == coef);
  }

  bool operator!=(const T& coef) const {
    return !(*this == coef);
  }

  //! Returns string representation using x_i as variables
  std::string toString() const;

  //! Returns string representation using var_i as variables
  std::string toString(const std::string& var) const;

  //! Returns string representation using vars as variables
  std::string toString(const std::vector<std::string>& vars) const;

  //! Returns string representation using vars as variables
  std::string toString(std::initializer_list<std::string> vars) const {
    return toString(std::vector<std::string>(vars));
  }

private:
  size_t dimension_;
  MonomialMap monomials_;
};

template <typename T, typename ExponentType, typename TermCompare>
MonomialSum<T, ExponentType, TermCompare>
operator+(MonomialSum<T, ExponentType, TermCompare> lhs, const MonomialSum<T, ExponentType, TermCompare>& rhs) {
  lhs += rhs;
  return lhs;
}

template <typename T, typename ExponentType, typename TermCompare>
MonomialSum<T, ExponentType, TermCompare>
operator+(MonomialSum<T, ExponentType, TermCompare> lhs, const Monomial<T, ExponentType>& monomial) {
  lhs += monomial;
  return lhs;
}

template <typename T, typename ExponentType, typename TermCompare>
MonomialSum<T, ExponentType, TermCompare>
operator+(const Monomial<T, ExponentType>& monomial, MonomialSum<T, ExponentType, TermCompare> rhs) {
  rhs += monomial;
  return rhs;
}

template <typename T, typename ExponentType, typename TermCompare>
MonomialSum<T, ExponentType, TermCompare>
operator+(MonomialSum<T, ExponentType, TermCompare> lhs, const Term<ExponentType>& term) {
  lhs += term;
  return lhs;
}

template <typename T, typename ExponentType, typename TermCompare>
MonomialSum<T, ExponentType, TermCompare>
operator+(const Term<ExponentType>& term, MonomialSum<T, ExponentType, TermCompare> rhs) {
  rhs += term;
  return rhs;
}

template <typename T, typename ExponentType, typename TermCompare>
MonomialSum<T, ExponentType, TermCompare> operator+(MonomialSum<T, ExponentType, TermCompare> lhs, const T& coef) {
  lhs += coef;
  return lhs;
}

template <typename T, typename ExponentType, typename TermCompare>
MonomialSum<T, ExponentType, TermCompare> operator+(const T& coef, MonomialSum<T, ExponentType, TermCompare> rhs) {
  rhs += coef;
  return rhs;
}

template <typename T, typename ExponentType, typename TermCompare>
MonomialSum<T, ExponentType, TermCompare>
operator-(MonomialSum<T, ExponentType, TermCompare> lhs, const MonomialSum<T, ExponentType, TermCompare>& rhs) {
  lhs -= rhs;
  return lhs;
}

template <typename T, typename ExponentType, typename TermCompare>
MonomialSum<T, ExponentType, TermCompare>
operator-(MonomialSum<T, ExponentType, TermCompare> lhs, const Monomial<T, ExponentType>& monomial) {
  lhs -= monomial;
  return lhs;
}

template <typename T, typename ExponentType, typename TermCompare>
MonomialSum<T, ExponentType, TermCompare>
operator-(const Monomial<T, ExponentType>& monomial, MonomialSum<T, ExponentType, TermCompare> rhs) {
  rhs *= -T(1);
  rhs += monomial;
  return rhs;
}

template <typename T, typename ExponentType, typename TermCompare>
MonomialSum<T, ExponentType, TermCompare>
operator-(MonomialSum<T, ExponentType, TermCompare> lhs, const Term<ExponentType>& term) {
  lhs -= term;
  return lhs;
}

template <typename T, typename ExponentType, typename TermCompare>
MonomialSum<T, ExponentType, TermCompare>
operator-(const Term<ExponentType>& term, MonomialSum<T, ExponentType, TermCompare> rhs) {
  rhs *= -T(1);
  rhs += term;
  return rhs;
}

template <typename T, typename ExponentType, typename TermCompare>
MonomialSum<T, ExponentType, TermCompare> operator-(MonomialSum<T, ExponentType, TermCompare> lhs, const T& coef) {
  lhs -= coef;
  return lhs;
}

template <typename T, typename ExponentType, typename TermCompare>
MonomialSum<T, ExponentType, TermCompare> operator-(const T& coef, MonomialSum<T, ExponentType, TermCompare> rhs) {
  rhs *= -T(1);
  rhs += coef;
  return rhs;
}

template <typename T, typename ExponentType, typename TermCompare>
MonomialSum<T, ExponentType, TermCompare>
operator*(MonomialSum<T, ExponentType, TermCompare> lhs, const MonomialSum<T, ExponentType, TermCompare>& rhs) {
  lhs *= rhs;
  return lhs;
}

template <typename T, typename ExponentType, typename TermCompare>
MonomialSum<T, ExponentType, TermCompare>
operator*(MonomialSum<T, ExponentType, TermCompare> lhs, const Monomial<T, ExponentType>& monomial) {
  lhs *= monomial;
  return lhs;
}

template <typename T, typename ExponentType, typename TermCompare>
MonomialSum<T, ExponentType, TermCompare>
operator*(const Monomial<T, ExponentType>& monomial, MonomialSum<T, ExponentType, TermCompare> rhs) {
  rhs *= monomial;
  return rhs;
}

template <typename T, typename ExponentType, typename TermCompare>
MonomialSum<T, ExponentType, TermCompare>
operator*(MonomialSum<T, ExponentType, TermCompare> lhs, const Term<ExponentType>& term) {
  lhs *= term;
  return lhs;
}

template <typename T, typename ExponentType, typename TermCompare>
MonomialSum<T, ExponentType, TermCompare>
operator*(const Term<ExponentType>& term, MonomialSum<T, ExponentType, TermCompare> rhs) {
  rhs *= term;
  return rhs;
}

template <typename T, typename ExponentType, typename TermCompare>
MonomialSum<T, ExponentType, TermCompare> operator*(MonomialSum<T, ExponentType, TermCompare> lhs, const T& coef) {
  lhs *= coef;
  return lhs;
}

template <typename T, typename ExponentType, typename TermCompare>
MonomialSum<T, ExponentType, TermCompare> operator*(const T& coef, MonomialSum<T, ExponentType, TermCompare> rhs) {
  rhs *= coef;
  return rhs;
}

template <
    typename T, typename ExponentType, typename TermCompare,
    typename = typename std::enable_if<std::is_signed<ExponentType>::value>::type>
MonomialSum<T, ExponentType, TermCompare>
operator/(MonomialSum<T, ExponentType, TermCompare> lhs, const Term<ExponentType>& term) {
  lhs /= term;
  return lhs;
}

template <typename T, typename ExponentType, typename TermCompare>
bool operator==(const Monomial<T, ExponentType>& monomial, const MonomialSum<T, ExponentType, TermCompare>& rhs) {
  return rhs == monomial;
}

template <typename T, typename ExponentType, typename TermCompare>
bool operator!=(const Monomial<T, ExponentType>& monomial, const MonomialSum<T, ExponentType, TermCompare>& rhs) {
  return !(rhs == monomial);
}

template <typename T, typename ExponentType, typename TermCompare>
bool operator==(const Term<ExponentType>& term, const MonomialSum<T, ExponentType, TermCompare>& rhs) {
  return rhs == term;
}

template <typename T, typename ExponentType, typename TermCompare>
bool operator!=(const Term<ExponentType>& term, const MonomialSum<T, ExponentType, TermCompare>& rhs) {
  return !(rhs == term);
}

template <typename T, typename ExponentType, typename TermCompare>
bool operator==(const T& coef, const MonomialSum<T, ExponentType, TermCompare>& rhs) {
  return rhs == coef;
}

template <typename T, typename ExponentType, typename TermCompare>
bool operator!=(const T& coef, const MonomialSum<T, ExponentType, TermCompare>& rhs) {
  return !(rhs == coef);
}

template <typename T, typename ExponentType, typename TermCompare>
void print(
    const MonomialSum<T, ExponentType, TermCompare>& p, const std::vector<std::string>& vars, std::ostream& out) {
  if (p.dimension() != vars.size()) {
    throw std::invalid_argument("Number of variables doesn't match dimension of the polynomial.");
  }

  if (p.isZero()) {
    out << "0";
    return;
  }

  bool is_first_term = true;

  for (const auto& it : p) {
    const auto& term = it.first;
    const auto& coef = it.second;

    if (!is_first_term) {
      out << " + ";
    }

    bool is_first_var = true;

    if (coef != T(1) || std::all_of(term.begin(), term.end(), [](const ExponentType& e) { return e == 0; })) {
      out << coef;
      is_first_var = false;
    }

    for (size_t i = 0; i < p.dimension(); ++i) {
      if (term[i] == 0) {
        continue;
      }

      if (!is_first_var) {
        out << " * ";
      }

      out << vars[i];

      if (term[i] != 1) {
        out << "^" << term[i];
      }

      is_first_var = false;
    }

    is_first_term = false;
  }
}

template <typename T, typename ExponentType, typename TermCompare>
void print(const MonomialSum<T, ExponentType, TermCompare>& p, const std::string& var, std::ostream& out) {
  std::vector<std::string> vars;
  vars.reserve(p.dimension());

  for (size_t i = 0; i < p.dimension(); ++i) {
    vars.push_back(var + "_" + std::to_string(i + 1));
  }

  print(p, vars, out);
}

template <typename T, typename ExponentType, typename TermCompare>
std::ostream& operator<<(std::ostream& out, const MonomialSum<T, ExponentType, TermCompare>& p) {
  print(p, "x", out);
  return out;
}

template <typename T, typename ExponentType, typename TermCompare>
std::string MonomialSum<T, ExponentType, TermCompare>::toString() const {
  return toString("x");
}

template <typename T, typename ExponentType, typename TermCompare>
std::string MonomialSum<T, ExponentType, TermCompare>::toString(const std::string& var) const {
  std::stringstream out;
  print(*this, var, out);
  return out.str();
}

template <typename T, typename ExponentType, typename TermCompare>
std::string MonomialSum<T, ExponentType, TermCompare>::toString(const std::vector<std::string>& vars) const {
  std::stringstream out;
  print(*this, vars, out);
  return out.str();
}

//! Evaluates a polynomial at a point.
template <typename T, typename ExponentType, typename TermCompare>
T evaluate(const MonomialSum<T, ExponentType, TermCompare>& m, const std::vector<T>& values);

//! Implementation for finite fields.
template <typename Ideal, typename ExponentType, typename TermCompare>
FieldElement<Ideal> evaluate(
    const MonomialSum<FieldElement<Ideal>, ExponentType, TermCompare>& m,
    const std::vector<FieldElement<Ideal>>& values) {
  if (m.dimension() != values.size()) {
    throw std::invalid_argument("The number of values provided doesn't match dimension of the polynomial.");
  }

  const auto zero = FieldElement<Ideal>(0);

  if (m.isZero()) {
    return zero;
  }

  FieldElement<Ideal> result(zero);

  for (const auto& it : m) {
    result += evaluate(it.second, it.first, values);
  }

  return result;
}

//! Standard polynomial uses unsigned integer type for variables' exponents
template <typename T>
using Polynomial = MonomialSum<T, uint32_t>;
using IntPolynomial = Polynomial<int>;

//! Laurent polynomial uses signed integer type for variables' exponents
template <typename T>
using LaurentPolynomial = MonomialSum<T, int32_t>;
using IntLaurentPolynomial = LaurentPolynomial<int>;
} // namespace polynomials
} // namespace crag

#endif // CRAG_POLYNOMIAL_H
