#include <gtest/gtest.h>

#include "fast_conjugacy_check.h"

#include "braid_group.h"
#include "random_word.h"

namespace crag {
namespace braidgroup {
namespace {

TEST(FastConjCheck, Ex_01) {
  using FF = finitefield::ZZ<10007>;

  const auto n = 16;

  FastConjugacyChecker<FF> checker(n);

  EXPECT_TRUE(checker.areNotConjugate("x1 x2 x3"_w, "x4"_w));
  EXPECT_FALSE(checker.areNotConjugate("x1 x2 x1"_w, "x2 x1 x2"_w)); 
  EXPECT_FALSE(checker.areNotConjugate("x4^-1 x5^-1 x2^-1 x1 x2 x1 x2 x5 x4"_w, "x2 x1 x2"_w)); 
  EXPECT_TRUE(checker.areNotConjugate("x4^-1 x5^-1 x2^-1 x1 x3 x1 x2 x5 x4"_w, "x2 x1 x2"_w)); 
}


TEST(FastConjCheck, Ex_02) {
  using FF = finitefield::ZZ<10007>;

  const auto n = 16;

  FastConjugacyChecker<FF> checker(n);

  std::mt19937_64 g(0);
  const size_t test_cases = 100;

  for (size_t i = 0; i < test_cases; ++i) {
    const auto random_u = random::randomWord(n - 1, 10, 20, g);
    const auto random_v = random::randomWord(n - 1, 10, 20, g);

    EXPECT_FALSE(checker.areNotConjugate(random_u, -random_v * random_u *random_v));
    EXPECT_TRUE(checker.areNotConjugate(random_u, random_v));
  }
}


} // namespace
} // namespace braidgroup
} // namespace crag
