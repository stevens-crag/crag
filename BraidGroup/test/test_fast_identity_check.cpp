#include <gtest/gtest.h>

#include "fast_identity_check.h"

#include "braid_group.h"
#include "random_word.h"

namespace crag {
namespace braidgroup {
namespace {

TEST(FastIdCheck, Ex_01) {
  using FF = finitefield::ZZ<199>;

  const auto n = 16;

  FastIdentityChecker<FF> id_checker(n);

  std::vector<Word> relations = getRelations(n);

  for (const auto& r : relations) {
    EXPECT_FALSE(id_checker.isNonTrivial(r));
  }
}

TEST(FastIdCheck, Ex_02) {
  using FF = finitefield::ZZ<199>;

  const auto n = 10;

  FastIdentityChecker<FF> id_checker(n);

  std::mt19937_64 g(0);
  const size_t test_cases = 100;

  for (size_t i = 0; i < test_cases; ++i) {
    const auto random_w = random::randomWord(n - 1, 10, 20, g);

    EXPECT_TRUE(id_checker.isNonTrivial(random_w));
  }
}

TEST(FastIdCheck, Hash) {
  using FF = finitefield::ZZ<199>;

  const size_t n = 16;
  BraidHasher<FF> hasher(n);

  std::mt19937 g(0);

  const auto random_w = random::randomWord(n - 1, 20, 30, g);
  const auto h = hasher(random_w);

  for (const auto& r : braidgroup::getRelations(n)) {
    EXPECT_EQ(h, hasher(random_w * r));
  }

  const auto other_random_w = random::randomWord(n - 1, 20, 30, g);

  EXPECT_NE(h, hasher(other_random_w));
}

TEST(FastIdCheck, HashVector) {
  using FF = finitefield::ZZ<199>;

  const size_t n = 16;

  BraidHasher<FF> hasher(n);
  std::mt19937 g(0);

  const std::vector<Word> random_words = {
      random::randomWord(n - 1, 20, 30, g),
      random::randomWord(n - 1, 20, 30, g),
      random::randomWord(n - 1, 20, 30, g),
  };

  const auto h = hasher(random_words);

  for (const auto& r : braidgroup::getRelations(n)) {
    std::vector<Word> other_words = {
        random_words[0] * r,
        random_words[1] * r,
        random_words[2] * r,
    };

    EXPECT_EQ(h, hasher(other_words));
  }
}
} // namespace
} // namespace braidgroup
} // namespace crag