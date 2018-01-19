#include "gtest/gtest.h"

#include "Word.h"

namespace crag {
namespace {

TEST(TestWord, TestInput1) {
  const auto w = "x1 x2^-2"_w;
  const Word result(std::vector<int>({1, -2, -2}));

  EXPECT_EQ(result, w);
}

} // namespace
} // namespace crag
