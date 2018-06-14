#include <gtest/gtest.h>
#include <random>

#include "LinkedBraidStructure.h"
#include "ShortBraidForm.h"
#include "random_word.h"

namespace crag {
namespace {

TEST(LinkedBraidStructure, Determinism) {
  const size_t n = 16;

  std::mt19937_64 g(0);

  const auto w = random::randomWord(n - 1, 100, 150, g);
  const auto w_short = shortBraidForm(n, w);

  for (size_t i = 0; i < 100; ++i) {
    EXPECT_EQ(w_short, shortBraidForm(n, w));
  }
}

TEST(LinkedBraidStructure, Equivalence) {
  const size_t n = 16;

  std::mt19937_64 g(0);

  for (size_t i = 0; i < 20; ++i) {
    const auto w = random::randomWord(n - 1, 200, 1000, g);
    const auto w_short = shortBraidForm(n, w);

    EXPECT_TRUE(areEqualBraids(n, w, w_short));
  }
}

TEST(LinkedBraidStructure, LongWords) {
  const size_t n = 16;

  std::mt19937_64 g(0);

  for (size_t i = 0; i < 5; ++i) {
    const auto w = random::randomWord(n - 1, 20000, g);
    const auto w_short = shortenBraid2(n, w);

    EXPECT_TRUE(areEqualBraids(n, w, w_short));
  }
}
} // namespace
} // namespace crag
