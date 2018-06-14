#include <gtest/gtest.h>

#include "colored_burau.h"

#include "braid_group.h"
#include "random_word.h"

namespace crag {
namespace coloredburau {
namespace {

using ZZ5 = finitefield::ZZ<5>;
using GF256 =
    finitefield::FieldElement<finitefield::IdealGeneratedByPolynomial<finitefield::ZZ<2>, 1, 1, 0, 1, 1, 0, 0, 0, 1>>;

TEST(ColoredBurau, B_3_matrix_mapping) {
  const LaurentPolynomial<ZZ5> zero(3);
  const auto unit = zero + ZZ5(1);

  const auto t_1 = zero + std::vector<int>{1, 0, 0};
  const auto t_2 = zero + std::vector<int>{0, 1, 0};
  const auto t_3 = zero + std::vector<int>{0, 0, 1};

  const auto t_2_inv = zero + std::vector<int>{0, -1, 0};
  const auto t_3_inv = zero + std::vector<int>{0, 0, -1};

  // clang-format off
  const auto b_1 = CBMatrix<ZZ5>(3, 1);
  EXPECT_EQ(Matrix<LaurentPolynomial<ZZ5>>(3, {
      -t_1, unit, zero,
      zero, unit, zero,
      zero, zero, unit,
  }), b_1);

  const auto b_1_inv = CBMatrix<ZZ5>(3, -1);
  EXPECT_EQ(Matrix<LaurentPolynomial<ZZ5>>(3, {
      -t_2_inv, t_2_inv, zero,
      zero, unit, zero,
      zero, zero, unit,
  }), b_1_inv);

  const auto b_2 = CBMatrix<ZZ5>(3, 2);
  EXPECT_EQ(Matrix<LaurentPolynomial<ZZ5>>(3, {
      unit, zero, zero,
      t_2, -t_2, unit,
      zero, zero, unit,
  }), b_2);

  const auto b_2_inv = CBMatrix<ZZ5>(3, -2);
  EXPECT_EQ(Matrix<LaurentPolynomial<ZZ5>>(3, {
      unit, zero, zero,
      unit, -t_3_inv, t_3_inv,
      zero, zero, unit,
  }), b_2_inv);
  // clang-format on
}

TEST(ColoredBurau, InducedPermutation) {
  EXPECT_EQ(Permutation({1, 0, 2}), permutation(3, Word(1)));
  EXPECT_EQ(Permutation({1, 0, 2}), permutation(3, Word(-1)));

  EXPECT_EQ(Permutation({0, 2, 1}), permutation(3, Word(2)));
  EXPECT_EQ(Permutation({0, 2, 1}), permutation(3, Word(-2)));

  EXPECT_EQ(Permutation(5), permutation(5, Word({1, 2, 3, 4, 4, -3, -2, -1})));

  EXPECT_EQ(Permutation(16), permutation(16, braidgroup::getPureBraidSubgroupGenerator(16, 1, 7)));
  EXPECT_EQ(Permutation(16), permutation(16, braidgroup::getPureBraidSubgroupGenerator(16, 3, 16)));
  EXPECT_EQ(Permutation(16), permutation(16, braidgroup::getPureBraidSubgroupGenerator(16, 5, 9)));
}

TEST(Permutations, Monomial) {
  const auto m = polynomials::IntStandardMonomial({1, 2, 3});

  // x y^2 z^3 goes to z x^2 y^3, and we store it like {2, 3, 1}
  EXPECT_EQ(polynomials::IntStandardMonomial({2, 3, 1}), permute(m, Permutation({2, 0, 1})));
}

TEST(Permutations, Polynomial) {
  const auto zero = LaurentPolynomial<ZZ5>(3);
  EXPECT_EQ(zero, permute(zero, Permutation({2, 1, 0})));

  const auto q = zero + polynomials::LaurentMonomial<ZZ5>({1, 2, 3}) + polynomials::LaurentMonomial<ZZ5>({-1, -2, -3});

  EXPECT_EQ(
      zero + polynomials::LaurentMonomial<ZZ5>({2, 3, 1}) + polynomials::LaurentMonomial<ZZ5>({-2, -3, -1}),
      permute(q, Permutation({2, 0, 1})));
}

TEST(Permutations, Matrix) {
  const LaurentPolynomial<ZZ5> zero(3);
  const auto unit = zero + ZZ5(1);

  const auto t_1 = zero + std::vector<int>{1, 0, 0};
  const auto t_2 = zero + std::vector<int>{0, 1, 0};
  const auto t_3 = zero + std::vector<int>{0, 0, 1};

  // clang-format off
  const auto b_1 = CBMatrix<ZZ5>(3, 1);
  EXPECT_EQ(Matrix<LaurentPolynomial<ZZ5>>(3, {
      -t_1, unit, zero,
      zero, unit, zero,
      zero, zero, unit,
  }), b_1);

  EXPECT_EQ(Matrix<LaurentPolynomial<ZZ5>>(3, {
      -t_3, unit, zero,
      zero, unit, zero,
      zero, zero, unit,
  }), permute(b_1, Permutation({2, 1, 0})));
  // clang-format on
}

TEST(ColoredBurau, B_3_mapping) {
  const auto b_1 = CBImage<ZZ5>(3, 1);
  const auto b_1_inv = CBImage<ZZ5>(3, -1);

  const auto b_2 = CBImage<ZZ5>(3, 2);
  const auto b_2_inv = CBImage<ZZ5>(3, -2);

  EXPECT_EQ(b_1 * b_2 * b_1, b_2 * b_1 * b_2);
  EXPECT_EQ(b_1 * b_2 * b_1, CBImage<ZZ5>(3, Word({1, 2, 1})));
  EXPECT_EQ(b_2 * b_1 * b_2, CBImage<ZZ5>(3, Word({2, 1, 2})));

  EXPECT_EQ(CBImage<ZZ5>(3, Word()), CBImage<ZZ5>(3, Word({1, 2, 1, -2, -1, -2})));
  EXPECT_EQ(CBImage<ZZ5>(3, Word()), b_1 * b_1_inv);
}

TEST(ColoredBurau, B_5_mapping) {
  const auto b_1 = CBImage<ZZ5>(5, 1);
  const auto b_2 = CBImage<ZZ5>(5, 2);
  const auto b_3 = CBImage<ZZ5>(5, 3);
  const auto b_4 = CBImage<ZZ5>(5, 4);

  EXPECT_EQ(b_1 * b_2 * b_1, b_2 * b_1 * b_2);
  EXPECT_EQ(b_2 * b_3 * b_2, b_3 * b_2 * b_3);
  EXPECT_EQ(b_3 * b_4 * b_3, b_4 * b_3 * b_4);

  EXPECT_EQ(b_1 * b_3, b_3 * b_1);
  EXPECT_EQ(b_1 * b_4, b_4 * b_1);
  EXPECT_EQ(b_2 * b_4, b_4 * b_2);
}

template <typename T>
void testBnMapping(size_t n) {
  for (int i = 1; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      if (j - i == 1) {
        EXPECT_EQ(CBImage<T>(n, Word({i, j, i})), CBImage<T>(n, Word({j, i, j})));
        EXPECT_EQ(CBImage<T>(n, Word()), CBImage<T>(n, Word({i, j, i, -j, -i, -j})));
      } else {
        EXPECT_EQ(CBImage<T>(n, Word({i, j})), CBImage<T>(n, Word({j, i})));
        EXPECT_EQ(CBImage<T>(n, Word()), CBImage<T>(n, Word({i, j, -i, -j})));
      }
    }
  }
}

TEST(ColoredBurau, B_n_mapping_over_zz5) {
  testBnMapping<ZZ5>(10);
}

TEST(ColoredBurau, B_n_mapping_over_gf256) {
  testBnMapping<GF256>(10);
}

TEST(ColoredBurau, EMultiplication1) {
  const auto m = Matrix<ZZ5>(3, {ZZ5(3), ZZ5(2), ZZ5(1), ZZ5(1), ZZ5(4), ZZ5(2), ZZ5(3), ZZ5(3), ZZ5(0)});
  const auto p = Permutation({2, 0, 1});
  const std::vector<ZZ5> t_values = {ZZ5(1), ZZ5(2), ZZ5(3)};

  const auto result = CBProjectionElement<ZZ5>(t_values, m, p) * Word(1);

  const auto expected = Matrix<ZZ5>(3, {ZZ5(1), ZZ5(0), ZZ5(1), ZZ5(2), ZZ5(0), ZZ5(2), ZZ5(1), ZZ5(1), ZZ5(0)});

  EXPECT_EQ(expected, result.matrix());
  EXPECT_EQ(Permutation({0, 2, 1}), result.permutation());
}

TEST(ColoredBurau, EMultiplication2) {
  const std::vector<ZZ5> t_values = {ZZ5(1), ZZ5(2), ZZ5(3)};

  const auto el = CBProjectionElement<ZZ5>(t_values);

  EXPECT_EQ(el, el * Word({1, 2, 1, -2, -1, -2}));
}

TEST(ColoredBurau, EMultiplication3) {
  using FF = GF256;

  const size_t n = 16;
  std::vector<FF> t_values(n, FF(0));

  std::mt19937 g(0);

  for (size_t i = 0; i < n; ++i) {
    t_values[i] = FF::random(g);
  }

  const auto test_cases = 100;

  for (size_t i = 0; i < test_cases; ++i) {
    const auto random_w = random::randomWord(n - 1, 20, 30, g);
    const auto cb_image = CBImage<FF>(n, random_w);

    EXPECT_EQ(project(cb_image, t_values), project(random_w, t_values)) << "case: i = " << i;
  }
}

TEST(ColoredBurau, OptimizedImageComputation) {
  using FF = GF256;

  const size_t n = 16;

  std::mt19937 g(0);

  for (size_t i = 0; i < 20; ++i) {
    const auto random_w = random::randomWord(n - 1, 20, 30, g);
    EXPECT_EQ(CBImageSlow<FF>(n, random_w), CBImage<FF>(n, random_w));
  }
}

TEST(ColoredBurau, Hash) {
  using FF = GF256;

  const size_t n = 16;
  std::vector<FF> t_values(n, FF(0));

  std::mt19937 g(0);

  for (size_t i = 0; i < n; ++i) {
    t_values[i] = FF::random(g);
  }

  const auto random_w = random::randomWord(n - 1, 20, 30, g);

  const auto braid_relations = braidgroup::getRelations(n);

  for (const auto& r : braid_relations) {
    EXPECT_EQ(
        std::hash<CBProjectionElement<FF>>()(project(random_w, t_values)),
        std::hash<CBProjectionElement<FF>>()(project(random_w * r, t_values)));
  }
}

TEST(Permutation, Ex_01) {
  const size_t n = 16;
  std::mt19937_64 g(0);

  const auto w1 = random::randomWord(n - 1, 20, 30, g);
  const auto w2 = random::randomWord(n - 1, 20, 30, g);

  const auto sigma1 = permutation(n, w1);
  const auto sigma2 = permutation(n, w2);

  EXPECT_EQ(permutation(n, w1 * w2), sigma1 * sigma2);

  EXPECT_EQ(Permutation(n), permutation(n, w1.inverse()) * sigma1);
  EXPECT_EQ(Permutation(n), sigma1 * permutation(n, w1.inverse()));
}
} // namespace
} // namespace coloredburau
} // namespace crag
