#include <gtest/gtest.h>
#include <random>

#include "random_subset.h"

namespace crag {
namespace random {
namespace {

TEST(RandomSubset, Ex_01) {
  std::mt19937_64 g(0);

  EXPECT_THROW({ subset(10, 0, 4, g); }, std::invalid_argument);
  EXPECT_THROW({ subset(1, 4, 2, g); }, std::invalid_argument);
}

TEST(RandomSubset, Ex_02) {
  std::mt19937_64 g(0);

  const auto result = subset(2, 0, 15, g);

  EXPECT_EQ(std::set<int>({2, 15}), result);
}

TEST(RandomSubset, Ex_03) {
  std::mt19937_64 g(0);

  const auto result = subsetVector(4, 1, 7, g);

  EXPECT_EQ(std::vector<int>({7, 2, 1, 5}), result);
}
} // namespace
} // namespace random
} // namespace crag