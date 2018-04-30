#include "gtest/gtest.h"

#include "Permutation.h"

namespace {

TEST(Permutation, DefaultConstructor) {
  Permutation p;

  EXPECT_EQ(0, p.size());
}

TEST(Permutation, Constructor) {
  Permutation p(2);

  ASSERT_EQ(2, p.size());

  EXPECT_EQ(0, p[0]);
  EXPECT_EQ(1, p[1]);

  EXPECT_EQ(Permutation({0, 1}), p);
}

TEST(Permutation, Print) {
  EXPECT_EQ("{}", Permutation().toString());
  EXPECT_EQ("{0, 1}", Permutation(2).toString());
  EXPECT_EQ("{1, 3, 2, 0}", Permutation({1, 3, 2, 0}).toString());
}

TEST(Permutation, Validation) {
  EXPECT_THROW({ Permutation(std::vector<int>{1, 1, 1}); }, std::invalid_argument);
  EXPECT_THROW({ Permutation({0, 1, 3}); }, std::invalid_argument);
}

TEST(Permutation, Multiplication) {
  EXPECT_EQ(Permutation(), Permutation() * Permutation());
  EXPECT_EQ(Permutation({}), Permutation({}) * Permutation({}));
  EXPECT_EQ(Permutation(5), Permutation() * Permutation(5));
  EXPECT_EQ(Permutation(2), Permutation({1, 0}) * Permutation({1, 0}));
  EXPECT_EQ(Permutation(5), Permutation(2) * Permutation(5));
  EXPECT_EQ(Permutation(7), Permutation(7) * Permutation(3));
  EXPECT_EQ(Permutation({0, 1, 3, 2}), Permutation({1, 0, 3, 2}) * Permutation({1, 0}));
  EXPECT_EQ(Permutation({0, 1, 4, 3, 2}), Permutation({1, 0}) * Permutation({1, 0, 4, 3, 2}));
}

TEST(Permutation, Inverse) {
  EXPECT_EQ(Permutation(), Permutation().inverse());
  EXPECT_EQ(Permutation({1, 0}), Permutation({1, 0}).inverse());
  EXPECT_EQ(Permutation({2, 0, 1}), Permutation({1, 2, 0}).inverse());
  EXPECT_EQ(Permutation({2, 0, 3, 1}), -Permutation({1, 3, 0, 2}));

  EXPECT_EQ(Permutation({1, 3, 0, 2}), -Permutation({1, 3, 0, 2}).inverse());
}

TEST(Permutation, Power) {
  EXPECT_EQ(Permutation(), Permutation().power(-3));
  EXPECT_EQ(Permutation(), Permutation().power(2));
  EXPECT_EQ(Permutation(2), Permutation({1, 0}).power(2));
  EXPECT_EQ(Permutation({2, 0, 1}), Permutation({1, 2, 0}).power(-1));
  EXPECT_EQ(Permutation(3), Permutation({1, 2, 0}).power(3));
  EXPECT_EQ(Permutation({3, 2, 1, 0}), Permutation({1, 3, 0, 2}).power(-2));
  EXPECT_EQ(Permutation(4), Permutation({1, 3, 0, 2}).power(0));
  EXPECT_EQ(Permutation({1, 3, 0, 2}), Permutation({1, 3, 0, 2}).power(1));
}

TEST(Permutation, IncreaseSize) {
  EXPECT_EQ(Permutation({1, 0}), Permutation({1, 0}).increaseSize(0));
  EXPECT_EQ(Permutation({1, 0, 2, 3, 4}), Permutation({1, 0}).increaseSize(5));
}

TEST(Permutation, GetHalfTwist) {
  EXPECT_EQ(Permutation(), Permutation::getHalfTwistPermutation(0));
  EXPECT_EQ(Permutation({3, 2, 1, 0}), Permutation::getHalfTwistPermutation(4));
}

TEST(Permutation, GetCycle) {
  EXPECT_EQ(Permutation(), Permutation::getCyclePermutation(0));
  EXPECT_EQ(Permutation({1, 2, 3, 0}), Permutation::getCyclePermutation(4));
}

TEST(Permutation, IsTrivial) {
  EXPECT_TRUE(Permutation().isTrivial());
  EXPECT_TRUE(Permutation(5).isTrivial());
  EXPECT_FALSE(Permutation({1, 3, 2, 0}).isTrivial());
}

TEST(Permutation, Geodesic) {
  EXPECT_EQ(std::vector<int>(), Permutation().geodesic());
  EXPECT_EQ(std::vector<int>(), Permutation(5).geodesic());
  EXPECT_EQ(std::vector<int>({0}), Permutation({1, 0}).geodesic());
  EXPECT_EQ(std::vector<int>({1, 0}), Permutation::getCyclePermutation(3).geodesic());
  EXPECT_EQ(std::vector<int>({2, 1, 0, 1, 2}), Permutation({3, 1, 2, 0}).geodesic());
}

TEST(Permutation, Difference) {
  EXPECT_EQ(std::numeric_limits<size_t>::max(), Permutation(2).difference(Permutation(3)));
  EXPECT_EQ(0, Permutation(3).difference(Permutation(3)));
  EXPECT_EQ(2, Permutation({1, 2, 0}).difference(Permutation({0, 2, 1})));
}

TEST(Permutation, ToCycles) {
  const std::vector<std::vector<int>> expected = {
      {0, 6, 4, 7},
      {1, 5, 3, 2, 9},
  };

  EXPECT_EQ(expected, toCycles(Permutation({6, 5, 9, 2, 7, 3, 4, 0, 8, 1})));
  EXPECT_TRUE(toCycles(Permutation(100)).empty());
}
} // namespace
