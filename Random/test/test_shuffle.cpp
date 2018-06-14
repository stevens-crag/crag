#include <gtest/gtest.h>
#include <random>

#include "shuffle.h"

namespace crag {
namespace random {
namespace {

TEST(Shuffle, Ex_01) {
  std::vector<int> items = {1, 2, 3, 4, 5};

  std::mt19937_64 g(0);
  random::shuffle(items.begin(), items.end(), g);

  EXPECT_EQ(std::vector<int>({3, 2, 5, 4, 1}), items);
}
} // namespace
} // namespace random
} // namespace crag