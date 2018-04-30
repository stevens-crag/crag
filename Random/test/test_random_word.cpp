#include <gtest/gtest.h>
#include <random>

#include "random_word.h"

namespace crag {
namespace random {
namespace {

TEST(RandomWord, Test1) {
  std::mt19937 g(1234);

  for (int i = 0; i < 10000; ++i) {
    auto w = randomWord(3, 100, g);

    ASSERT_EQ(100, w.length());
  }
}

TEST(RandomWord, Test2) {
  std::mt19937 g(1234);

  ASSERT_EQ(4, randomWord(1, 4, g).length());
}

TEST(RandomWord, Test3) {
  std::mt19937_64 g(0);

  EXPECT_EQ(Word(), randomWord(Permutation(10), g));
  EXPECT_EQ(Word(1), randomWord(Permutation({1, 0}), g));
  EXPECT_EQ(Word(-1), randomWord(Permutation({1, 0}), g));
}
} // namespace
} // namespace random
} // namespace crag
