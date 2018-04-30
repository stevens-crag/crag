#include <gtest/gtest.h>
#include <random>

#include "random_permutation.h"

namespace crag {
namespace random {
namespace {

TEST(RandomPermutation, Ex_01) {
  std::mt19937_64 g(0);
  EXPECT_EQ(Permutation({6, 5, 9, 2, 7, 3, 4, 0, 8, 1}), randomPermutation(10, g));
}
} // namespace
} // namespace random
} // namespace crag
