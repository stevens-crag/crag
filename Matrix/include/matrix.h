#pragma once

#ifndef CRAG_MATRIX_H
#define CRAG_MATRIX_H

#include <ostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

namespace crag {
namespace matrix {

//! Plain template matrix over (discrete) ring
template <typename T>
class Matrix {
public:
  Matrix() = delete;

  //! Creates (n, n) matrix filled with T(0)
  explicit Matrix(size_t n)
      : rows_(n, std::vector<T>(n, T(0))) {
    if (n == 0) {
      throw std::invalid_argument("Dimensions of the matrix must be non zero.");
    }
  }

  //! Creates (n, n) matrix filled with value
  Matrix(size_t n, const T& value)
      : rows_(n, std::vector<T>(n, value)) {
    if (n == 0) {
      throw std::invalid_argument("Dimensions of the matrix must be non zero.");
    }
  }

  //! Creates (n, n) matrix filled with values
  Matrix(size_t n, const std::vector<T>& values)
      : rows_(toRows_(n, n, values)) {
    if (n == 0) {
      throw std::invalid_argument("Dimensions of the matrix must be non zero.");
    }
  }

  //! Creates (n, m) matrix filled with T(0)
  explicit Matrix(std::pair<size_t, size_t> size)
      : rows_(size.first, std::vector<T>(size.second, T(0))) {
    if ((size.first == 0) || (size.second == 0)) {
      throw std::invalid_argument("Dimensions of the matrix must be non zero.");
    }
  }

  //! Creates (n, m) matrix filled with value
  Matrix(std::pair<size_t, size_t> size, const T& value)
      : rows_(size.first, std::vector<T>(size.second, value)) {
    if ((size.first == 0) || (size.second == 0)) {
      throw std::invalid_argument("Dimensions of the matrix must be non zero.");
    }
  }

  //! Creates (n, m) matrix filled with values
  Matrix(std::pair<size_t, size_t> size, const std::vector<T>& values)
      : rows_(toRows_(size.first, size.second, values)) {
    if ((size.first == 0) || (size.second == 0)) {
      throw std::invalid_argument("Dimensions of the matrix must be non zero.");
    }
  }

  //! Returns the number of rows
  size_t size1() const {
    return rows_.size();
  }

  //! Returns the number of columns
  size_t size2() const {
    return rows_.front().size();
  }

  //! Exchanges the content of matrices
  void swap(Matrix& other) {
    std::swap(rows_, other.rows_);
  }

  T& operator()(size_t i, size_t j) {
    return rows_[i][j];
  }

  const T& operator()(size_t i, size_t j) const {
    return rows_[i][j];
  }

  bool operator==(const Matrix& other) const {
    if ((size1() != other.size1()) || (size2() != other.size2())) {
      return false;
    }

    for (size_t i = 0; i < size1(); ++i) {
      for (size_t j = 0; j < size2(); ++j) {
        if (rows_[i][j] != other.rows_[i][j]) {
          return false;
        }
      }
    }

    return true;
  }

  bool operator!=(const Matrix& other) const {
    return !(*this == other);
  }

  bool operator==(const T& val) const {
    return (size1() == 1) && (size2() == 1) && (rows_[0][0] == val);
  }

  bool operator!=(const T& val) const {
    return !(*this == val);
  }

  Matrix& operator*=(const T& coef) {
    for (size_t i = 0; i < size1(); ++i) {
      for (size_t j = 0; j < size2(); ++j) {
        rows_[i][j] *= coef;
      }
    }

    return *this;
  }

  Matrix& operator*=(const Matrix& other) {
    if ((size1() != other.size1()) || (size2() != other.size2())) {
      throw std::invalid_argument("Dimensions of matrices don't match.");
    }

    std::vector<std::vector<T>> new_rows;
    new_rows.reserve(size1());

    for (size_t i = 0; i < size1(); ++i) {
      std::vector<T> row;
      row.reserve(other.size2());

      for (size_t j = 0; j < other.size2(); ++j) {
        T m_ij = rows_[i][0] * other.rows_[0][j];

        for (size_t k = 1; k < size2(); ++k) {
          m_ij += rows_[i][k] * other.rows_[k][j];
        }

        row.push_back(std::move(m_ij));
      }

      new_rows.push_back(std::move(row));
    }

    rows_ = std::move(new_rows);

    return *this;
  }

  Matrix& operator+=(const Matrix& other) {
    if ((size1() != other.size1()) || (size2() != other.size2())) {
      throw std::invalid_argument("Dimensions of matrices don't match.");
    }

    for (size_t i = 0; i < size1(); ++i) {
      for (size_t j = 0; j < size2(); ++j) {
        rows_[i][j] += other.rows_[i][j];
      }
    }

    return *this;
  }

  Matrix& operator-=(const Matrix& other) {
    if ((size1() != other.size1()) || (size2() != other.size2())) {
      throw std::invalid_argument("Dimensions of matrices don't match.");
    }

    for (size_t i = 0; i < size1(); ++i) {
      for (size_t j = 0; j < size2(); ++j) {
        rows_[i][j] -= other.rows_[i][j];
      }
    }

    return *this;
  }

  std::string toString() const;

private:
  std::vector<std::vector<T>> rows_;

  std::vector<std::vector<T>> toRows_(size_t n, size_t m, const std::vector<T>& values) {
    if ((n * m) != values.size()) {
      throw std::invalid_argument("Dimensions of the matrix doesn't match the number of values provided.");
    }

    std::vector<std::vector<T>> rows;
    rows.reserve(n);

    auto it = values.begin();

    for (size_t i = 0; i < n; ++i) {
      rows.push_back(std::vector<T>(it, it + m));
      it += m;
    }

    return rows;
  }
};

template <typename T>
bool operator==(const T& val, const Matrix<T>& m) {
  return m == val;
}

template <typename T>
bool operator!=(const T& val, const Matrix<T>& m) {
  return m != val;
}

template <typename T>
Matrix<T> operator*(Matrix<T> lhs, const T& coef) {
  return lhs *= coef;
}

template <typename T>
Matrix<T> operator*(const T& coef, Matrix<T> m) {
  return m *= coef;
}

template <typename T>
Matrix<T> operator*(Matrix<T> lhs, const Matrix<T>& rhs) {
  return lhs *= rhs;
}

template <typename T>
Matrix<T> operator+(Matrix<T> lhs, const Matrix<T>& rhs) {
  return lhs += rhs;
}

template <typename T>
Matrix<T> operator-(Matrix<T> lhs, const Matrix<T>& rhs) {
  return lhs -= rhs;
}

// boost::numeric::ublas style matrix output
template <typename T>
std::ostream& operator<<(std::ostream& out, const Matrix<T>& m) {
  out << "[" << m.size1() << "," << m.size2() << "](";

  for (size_t i = 0; i < m.size1(); ++i) {
    if (i > 0) {
      out << ",";
    }

    out << "(";

    for (size_t j = 0; j < m.size2(); ++j) {
      if (j > 0) {
        out << ",";
      }

      out << m(i, j);
    }

    out << ")";
  }

  out << ")";

  return out;
}

template <typename T>
std::string Matrix<T>::toString() const {
  std::stringstream out;
  out << *this;
  return out.str();
}

template <typename T>
Matrix<T> zero(size_t n) {
  return Matrix<T>(n);
}

template <typename T>
Matrix<T> zero(size_t n, const T& zero) {
  return Matrix<T>(n, zero);
}

template <typename T>
Matrix<T> zero(size_t n, size_t m) {
  return Matrix<T>(n, m);
}

template <typename T>
Matrix<T> zero(size_t n, size_t m, const T& zero) {
  return Matrix<T>(n, m, zero);
}

template <typename T>
Matrix<T> eye(size_t n, const T& zero, const T& value) {
  Matrix<T> m(n, zero);

  for (size_t i = 0; i < n; ++i) {
    m(i, i) = value;
  }

  return m;
}

template <typename T>
Matrix<T> eye(size_t n, const T& value) {
  return eye<T>(n, T(0), value);
}

template <typename T>
Matrix<T> eye(size_t n, const T& zero, const std::vector<T>& values) {
  if (n != values.size()) {
    throw std::invalid_argument("Dimension of the matrix doesn't match the number of values provided.");
  }

  Matrix<T> m(n, zero);

  for (size_t i = 0; i < n; ++i) {
    m(i, i) = values[i];
  }

  return m;
}

template <typename T>
Matrix<T> eye(size_t n, const std::vector<T>& values) {
  return eye(n, T(0), values);
}

template <typename T>
Matrix<T> eye(const T& zero, const std::vector<T>& values) {
  return eye(values.size(), zero, values);
}

template <typename T>
Matrix<T> eye(const std::vector<T>& values) {
  return eye(values.size(), values);
}

template <typename T>
Matrix<T> unit(size_t n) {
  return eye<T>(n, T(0), T(1));
}

template <typename T>
Matrix<T> unit(size_t n, const T& zero, const T& unit) {
  return eye<T>(n, zero, unit);
}

template <typename T>
Matrix<T> matrix(size_t n, size_t m, const std::vector<T>& values) {
  return Matrix<T>(n, m, std::move(values));
}

template <typename T>
Matrix<T> matrix(size_t n, const std::vector<T>& values) {
  return Matrix<T>(n, std::move(values));
}

template <typename T>
bool isSquare(const Matrix<T>& m) {
  return m.size1() == m.size2();
}

template <typename T> T tr(const Matrix<T> &m) {
  if (m.size1() != m.size2()) {
    throw std::invalid_argument("Dimensions of the matrix don't match.");
  }

  T result(m(0, 0));

  for (size_t i = 1; i < m.size1(); ++i) {
    result += m(i, i);
  }

  return result;
}

template <typename T> std::vector<T> charPoly(const Matrix<T> &a) {
  if (a.size1() != a.size2()) {
    throw std::invalid_argument("Dimensions of the matrix don't match.");
  }

  size_t n = a.size1();

  std::vector<T> result;

  auto m = zero<T>(n);

  const auto id = unit<T>(n);

  result.push_back(T(1));

  for (size_t k = 1; k <= n; ++k) {
    m = a * (m + result.back() * id);
    result.push_back((T(-1) / T(k)) * tr(m));
  }

  std::reverse(result.begin(), result.end());

  return result;
}

} // namespace matrix
} // namespace crag

#endif // CRAG_MATRIX_H
