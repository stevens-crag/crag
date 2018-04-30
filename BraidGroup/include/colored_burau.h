#pragma once

#ifndef CRAG_COLORED_BURAU_H
#define CRAG_COLORED_BURAU_H

#include <ostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "FiniteField.h"
#include "Permutation.h"
#include "Word.h"
#include "matrix.h"
#include "polynomial.h"

namespace crag {
namespace coloredburau {

using finitefield::FieldElement;
using matrix::Matrix;
using polynomials::LaurentPolynomial;
using polynomials::Monomial;
using polynomials::MonomialSum;
using polynomials::Term;

//! Computes matrix part of the image of a generator x_i^{\pm 1} of B_n in the colored Burau group.
//! Uses type T as a coefficient ring for Laurent polynomials in n variables.
//! Requires -n < i < n.
template <typename T>
Matrix<LaurentPolynomial<T>> CBMatrix(size_t n, int i) {
  if ((i == 0) || (static_cast<size_t>(std::abs(i)) >= n)) {
    throw std::invalid_argument("Index of a generator is out of range.");
  }

  const LaurentPolynomial<T> zero(n);
  const auto unit = zero + T(1);

  auto m = matrix::unit<LaurentPolynomial<T>>(n, zero, unit);

  if (i > 0) {
    std::vector<typename LaurentPolynomial<T>::exponent_t> term_t_i(n, 0);
    term_t_i[i - 1] = 1;

    const auto t_i = zero + term_t_i;

    if (i > 1) {
      m(i - 1, i - 2) = t_i;
    }

    m(i - 1, i - 1) = -t_i;
    m(i - 1, i) = unit;
  } else {
    i = -i;

    std::vector<typename LaurentPolynomial<T>::exponent_t> term_t_i_inverse(n, 0);
    term_t_i_inverse[i] = -1;

    const auto t_i_inverse = zero + term_t_i_inverse;

    if (i > 1) {
      m(i - 1, i - 2) = unit;
    }

    m(i - 1, i - 1) = -t_i_inverse;
    m(i - 1, i) = t_i_inverse;
  }

  return m;
}

template <typename T, typename ExponentType>
Monomial<T, ExponentType> permute(const T& coef, const Term<ExponentType>& term, const Permutation& p) {
  if (term.size() != p.size()) {
    throw std::invalid_argument("Dimensions of the term and the permutation don't match.");
  }

  auto new_term = term;

  for (size_t i = 0; i < term.size(); ++i) {
    new_term[p[i]] = term[i];
  }

  return Monomial<T, ExponentType>(coef, std::move(new_term));
}

template <typename T, typename ExponentType>
Monomial<T, ExponentType> permute(const Monomial<T, ExponentType>& m, const Permutation& p) {
  return permute(m.coef(), m.term(), p);
}

template <typename T, typename ExponentType, typename TermCompare>
MonomialSum<T, ExponentType, TermCompare>
permute(const MonomialSum<T, ExponentType, TermCompare>& polynomial, const Permutation& p) {
  if (polynomial.dimension() != p.size()) {
    throw std::invalid_argument("Dimensions of the polynomial and the permutation don't match.");
  }

  if (polynomial.isZero()) {
    return polynomial;
  }

  auto it = polynomial.begin();

  MonomialSum<T, ExponentType, TermCompare> result(permute(it->second, it->first, p));

  while (++it != polynomial.end()) {
    result += permute(it->second, it->first, p);
  }

  return result;
}

template <typename T, typename ExponentType, typename TermCompare>
Matrix<MonomialSum<T, ExponentType, TermCompare>>
permute(Matrix<MonomialSum<T, ExponentType, TermCompare>> m, const Permutation& p) {
  for (size_t i = 0; i < m.size1(); ++i) {
    for (size_t j = 0; j < m.size2(); ++j) {
      m(i, j) = permute(m(i, j), p);
    }
  }

  return m;
}

//! Returns the permutation p \in S_n corresponding to a braid word w \in B_n
Permutation permutation(size_t n, const Word& w);

//! Represents element of colored Burau group
template <typename T>
class CBElement {
public:
  using matrix_t = Matrix<LaurentPolynomial<T>>;

  CBElement(matrix_t m, Permutation p)
      : matrix_(std::move(m))
      , permutation_(std::move(p)) {
    if (!matrix::isSquare(matrix_) || (matrix_.size1() != permutation_.size())) {
      throw std::invalid_argument("Invalid matrix/permutation pair for colored Burau element");
    }
  }

  bool operator==(const CBElement& other) const {
    return (matrix_ == other.matrix_) && (permutation_ == other.permutation_);
  }

  const matrix_t& matrix() const {
    return matrix_;
  }

  const Permutation& permutation() const {
    return permutation_;
  }

  //! n from B_n
  size_t n() const {
    return matrix_.size1();
  }

  CBElement& operator*=(const CBElement& other) {
    if (n() != other.n()) {
      throw std::invalid_argument("Dimensions don't match.");
    }

    matrix_ *= permute(other.matrix_, permutation_);

    // change order of permutations' product since out operator* is not the canonical product
    permutation_ = other.permutation_ * permutation_;

    return *this;
  }

  //! Optimized multiplication by a word w
  CBElement& operator*=(const Word& w) {
    for (auto i : w) {
      if ((i == 0) || (static_cast<size_t>(std::abs(i)) >= n())) {
        throw std::invalid_argument("Index of a generator is out of range.");
      }

      if (i > 0) {
        // permutation replaces t_i with t_{sigma(i)}
        typename LaurentPolynomial<T>::term_t t(n(), 0);
        t[permutation_[i - 1]] = 1;

        const typename LaurentPolynomial<T>::monomial_t m(-T(1), t);

        for (size_t j = 0; j < n(); ++j) {
          matrix_(j, i) += matrix_(j, i - 1);
          matrix_(j, i - 1) *= m;

          if (i > 1) {
            matrix_(j, i - 2) -= matrix_(j, i - 1);
          }
        }
      } else {
        i = -i;

        // permutation replaces t_{i+1} with t_{sigma(i+1)}
        typename LaurentPolynomial<T>::term_t t(n(), 0);
        t[permutation_[i]] = -1;

        const typename LaurentPolynomial<T>::monomial_t m(-T(1), t);

        for (size_t j = 0; j < n(); ++j) {
          if (i > 1) {
            matrix_(j, i - 2) += matrix_(j, i - 1);
          }

          matrix_(j, i - 1) *= m;
          matrix_(j, i) -= matrix_(j, i - 1);
        }
      }

      permutation_.change(i - 1, i);
    }

    return *this;
  }

  std::string toString() const;

private:
  matrix_t matrix_;
  Permutation permutation_;
};

template <typename T>
CBElement<T> operator*(CBElement<T> lhs, const CBElement<T>& rhs) {
  return lhs *= rhs;
}

template <typename T>
std::ostream& operator<<(std::ostream& out, const CBElement<T>& cb_element) {
  out << "(" << cb_element.matrix() << ", " << cb_element.permutation() << ")";
  return out;
}

template <typename T>
std::string CBElement<T>::toString() const {
  std::stringstream out;
  out << *this;
  return out.str();
}

//! Computes the image of a generator x_i^{\pm 1} of B_n in the colored Burau group.
//! Uses type T as a coefficient ring for Laurent polynomials in n variables.
//! Requires -n < i < n.
template <typename T>
CBElement<T> CBImage(size_t n, int i) {
  return CBElement<T>(CBMatrix<T>(n, i), permutation(n, Word(i)));
}

//! Computes the image of w \in B_n in the colored Burau group.
//! Uses type T as a coefficient ring for Laurent polynomials in n variables.
//! Use with caution: for large w the resulting matrix may contain very large polynomials.
template <typename T>
CBElement<T> CBImage(size_t n, const Word& w) {
  const LaurentPolynomial<T> zero(n);
  const auto unit = zero + T(1);

  auto result = CBElement<T>(matrix::unit<LaurentPolynomial<T>>(n, zero, unit), Permutation(n));

  return result *= w;
}

//! Slow (not optimized, but more clear) version
template <typename T>
CBElement<T> CBImageSlow(size_t n, const Word& w) {
  const LaurentPolynomial<T> zero(n);
  const auto unit = zero + T(1);

  auto result = CBElement<T>(matrix::unit<LaurentPolynomial<T>>(n, zero, unit), Permutation(n));

  for (const auto i : w) {
    result *= CBImage<T>(n, i);
  }

  return result;
}

template <typename T, typename ExponentType, typename TermCompare>
Matrix<T> evaluate(const Matrix<MonomialSum<T, ExponentType, TermCompare>>& m, const std::vector<T>& values);

template <typename Ideal, typename ExponentType, typename TermCompare>
Matrix<FieldElement<Ideal>> evaluate(
    const Matrix<MonomialSum<FieldElement<Ideal>, ExponentType, TermCompare>>& m,
    const std::vector<FieldElement<Ideal>>& values) {
  Matrix<FieldElement<Ideal>> result(std::make_pair(m.size1(), m.size2()));

  for (size_t i = 0; i < m.size1(); ++i) {
    for (size_t j = 0; j < m.size2(); ++j) {
      result(i, j) = polynomials::evaluate(m(i, j), values);
    }
  }

  return result;
}

//! Represents an element of a finite set of pairs (M, sigma) parametrized by t-values (tau_1,...,tau_n) \in T^n,
//! where M is an (n, n) matrix over ring T, and sigma is a permutation of size n.
//! Colored Burau group corresponding to B_n acts on this set, and this action is called E-multiplication.
template <typename T>
class CBProjectionElement {
public:
  using matrix_t = Matrix<T>;

  //! Constructs a unit, i.e. (unit matrix, unit permutation)
  explicit CBProjectionElement(std::vector<T> t_values)
      : t_values_(std::move(t_values))
      , matrix_(matrix::unit<T>(t_values_.size()))
      , permutation_(t_values_.size()) {
    if (t_values_.size() < 2) {
      throw std::invalid_argument("Require at least 2 t-values.");
    }
  }

  CBProjectionElement(std::vector<T> t_values, matrix_t m, Permutation p)
      : t_values_(std::move(t_values))
      , matrix_(std::move(m))
      , permutation_(std::move(p)) {
    const auto size = t_values_.size();

    if (size < 2) {
      throw std::invalid_argument("Require at least 2 t-values.");
    }

    if (!matrix::isSquare(matrix_) || (size != matrix_.size1()) || (size != permutation_.size())) {
      throw std::invalid_argument("Dimensions of t-values/matrix/permutation don't match.");
    }
  }

  const std::vector<T>& tValues() const {
    return t_values_;
  }

  const matrix_t& matrix() const {
    return matrix_;
  }

  const Permutation& permutation() const {
    return permutation_;
  }

  size_t n() const {
    return t_values_.size();
  }

  //! E-multiplication
  CBProjectionElement& operator*=(const CBElement<T>& cb_element) {
    if (n() != cb_element.n()) {
      throw std::invalid_argument("Dimensions of matrices don't match.");
    }

    matrix_ *= evaluate(permute(cb_element.matrix(), permutation_), t_values_);

    // change order of permutations' product since out operator* is not the canonical product
    permutation_ = cb_element.permutation() * permutation_;

    return *this;
  }

  //! E-multiplication by a braid word w, performs multiplication iteratively for each generator of w,
  //! so it doesn't compute the image of w in colored Burau group.
  //! Multiplication by each generator is optimized to modify only 3 columns of the original matrix.
  CBProjectionElement& operator*=(const Word& w) {
    const auto n = this->n();

    for (const auto i : w) {
      const auto index = std::abs(i);

      if ((index < 1) || (index + 1 > n)) {
        throw std::invalid_argument("Cannot perform E-multiplication by w, generator's index is out of range.");
      }

      if (i > 0) {
        const auto& t = t_values_[permutation_[index - 1]];

        for (size_t j = 0; j < n; ++j) {
          matrix_(j, index) += matrix_(j, index - 1);
          matrix_(j, index - 1) *= -t;

          if (index > 1) {
            matrix_(j, index - 2) -= matrix_(j, index - 1);
          }
        }
      } else {
        const auto& t = t_values_[permutation_[index]].inverse();

        for (size_t j = 0; j < n; ++j) {
          if (index > 1) {
            matrix_(j, index - 2) += matrix_(j, index - 1);
          }

          matrix_(j, index - 1) *= -t;
          matrix_(j, index) -= matrix_(j, index - 1);
        }
      }

      permutation_.change(index - 1, index);
    }

    return *this;
  }

  std::string toString() const;

  bool operator==(const CBProjectionElement& other) const {
    return (t_values_ == other.t_values_) && (matrix_ == other.matrix_) && (permutation_ == other.permutation_);
  }

  bool operator!=(const CBProjectionElement& other) const {
    return !(*this == other);
  }

private:
  std::vector<T> t_values_;
  matrix_t matrix_;
  Permutation permutation_;
};

template <typename T>
CBProjectionElement<T> operator*(CBProjectionElement<T> lhs, const CBElement<T>& cb_el) {
  return lhs *= cb_el;
}

template <typename T>
CBProjectionElement<T> operator*(CBProjectionElement<T> lhs, const Word& w) {
  return lhs *= w;
}

template <typename T>
std::ostream& operator<<(std::ostream& out, const CBProjectionElement<T>& el) {
  out << "[";

  bool is_first = true;
  for (const auto& t_val : el.tValues()) {
    if (!is_first) {
      out << ",";
    }

    out << t_val;

    is_first = false;
  }

  out << "]";

  out << "(" << el.matrix() << ", " << el.permutation() << ")";

  return out;
}

template <typename T>
std::string CBProjectionElement<T>::toString() const {
  std::stringstream out;
  out << *this;
  return out.str();
}

//! Acts on the trivial pair (E, id) by cb_el using provided t-values.
template <typename T>
CBProjectionElement<T> project(const CBElement<T>& cb_el, std::vector<T> t_values) {
  CBProjectionElement<T> result(std::move(t_values));
  return result *= cb_el;
}

//! Acts on the trivial pair (E, id) by w using provided t-values.
template <typename T>
CBProjectionElement<T> project(const Word& w, std::vector<T> t_values) {
  CBProjectionElement<T> result(std::move(t_values));
  return result *= w;
}
} // namespace coloredburau
} // namespace crag

namespace std {

template <typename T>
struct hash<crag::coloredburau::CBProjectionElement<T>> {
public:
  size_t operator()(const crag::coloredburau::CBProjectionElement<T>& element) const {
    return std::hash<std::string>()(element.toString());
  }
};
} // namespace std
#endif // CRAG_COLORED_BURAU_H
