#include <gtest/gtest.h>

#include "cloaking_element.h"

namespace crag {
namespace coloredburau {
namespace {

using finitefield::ZZ;
using ZZ5 = ZZ<5>;

TEST(CloakingElement, Ex_01) {
  const size_t n = 16;
  const size_t a = 2;
  const size_t b = 4;

  // clang-format off
  const std::vector<ZZ5> t_values = {
      ZZ5(3), ZZ5(2), ZZ5(1), ZZ5(4), ZZ5(1), ZZ5(2), ZZ5(3), ZZ5(2),
      ZZ5(3), ZZ5(2), ZZ5(4), ZZ5(4), ZZ5(2), ZZ5(2), ZZ5(3), ZZ5(3),
  };
  // clang-format on

  std::mt19937 g(0);

  const size_t test_cases = 50;

  for (size_t i = 0; i < test_cases; ++i) {
    const auto elem = random::randomWord(n - 1, 10, 20, g);
    const auto sigma = coloredburau::permutation(n, elem);

    // cloaking elements with x_i^2 in the middle and t_a = t_b = 1
    const auto c2 = generateCloakingElement(10, n, a, b, 20, 30, sigma, 2, g);

    // cloaking elements with x_i^4 in the middle and t_a, t_b such that t_a * t_b = -1
    const auto c4 = generateCloakingElement(10, n, 1, 5, 20, 30, sigma, 4, g);

    const auto expected = coloredburau::project(elem, t_values);

    EXPECT_EQ(Permutation(n), sigma * expected.permutation());

    EXPECT_EQ(expected, coloredburau::project(elem * c2, t_values));
    EXPECT_EQ(expected, coloredburau::project(elem * c4, t_values));
  }
}

TEST(CloakingElement, WordByPermutation) {
  const size_t n = 16;
  std::mt19937_64 g(0);

  for (size_t i = 0; i < 10; ++i) {
    const auto sigma = random::randomPermutation(n, g);
    const auto w = random::randomWord(sigma, g);

    EXPECT_EQ(sigma, coloredburau::permutation(n, w));
  }
}

TEST(CloakingElement, Canonical) {
  const size_t n = 16;
  const size_t a = 2;
  const size_t b = 4;
  const size_t L = 30;

  // clang-format off
  const std::vector<ZZ5> t_values = {
      ZZ5(3), ZZ5(2), ZZ5(1), ZZ5(4), ZZ5(1), ZZ5(2), ZZ5(3), ZZ5(2),
      ZZ5(3), ZZ5(2), ZZ5(4), ZZ5(4), ZZ5(2), ZZ5(2), ZZ5(3), ZZ5(3),
  };
  // clang-format on

  std::mt19937_64 g(0);

  const size_t test_cases = 50;

  for (size_t i = 0; i < test_cases; ++i) {
    const auto elem = random::randomWord(n - 1, 10, 20, g);
    const auto sigma = coloredburau::permutation(n, elem);

    // cloaking elements with x_i^2 in the middle and t_a = t_b = 1
    const auto c2 = generateCanonicalCloakingEl(n, a, b, L, sigma, 2, g);

    // cloaking elements with x_i^4 in the middle and t_a, t_b such that t_a * t_b = -1
    const auto c4 = generateCanonicalCloakingEl(n, a, b, L, sigma, 4, g);

    const auto expected = coloredburau::project(elem, t_values);

    EXPECT_EQ(Permutation(n), sigma * expected.permutation());

    EXPECT_EQ(expected, coloredburau::project(elem * c2, t_values));
    EXPECT_EQ(expected, coloredburau::project(elem * c4, t_values));
    EXPECT_EQ(expected, expected * c2);
    EXPECT_EQ(expected, expected * c4);
  }
}
} // namespace
} // namespace coloredburau
} // namespace crag
