#include <gtest/gtest.h>
#include <random>
#include <algorithm>
#include "parallel.h"

namespace crag {
namespace parallel {
namespace {

TEST(ParallelMap, Ex_01) {
  const size_t size = 1 << 28;

  std::vector<int> items(size);

  std::mt19937 g(0);

  for (size_t i = 0; i < size; ++i) {
    items[i] = g() % 1000;
  }

  const auto squares = map(items, [](int i) {
    return i * i;
  });

  for (size_t i = 0; i < size; ++i) {
    EXPECT_EQ(items[i] * items[i], squares[i]);
  }
}

TEST(ParallelMap, Ex_02) {
  const size_t size = 1 << 20;

  std::vector<int> items(size);

  std::mt19937 g(0);

  for (size_t i = 0; i < size; ++i) {
    items[i] = g() % 1000;
  }

  const auto bools = bmap(items, [](int i) {
    return i >= 0;
  });

  EXPECT_TRUE(std::all_of(bools.begin(), bools.end(), [] (bool b) {
    return b;
  }));
}
}
}
}
