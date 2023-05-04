#include <gtest/gtest.h>
#include <numeric>
#include <random>

#include "LinkedBraidStructure.h"

#include "random_word.h"
#include "stochastic_rewrite.h"

namespace crag {
namespace stochasticrewrite {
namespace {

TEST(StochasticRewriteTest, TestRandomPartition) {
  std::mt19937 g(1234);

  const auto result = randomPartition(16, 3, 10, g);

  EXPECT_EQ(16, std::accumulate(result.begin(), result.end(), 0));
  EXPECT_TRUE(std::all_of(result.begin(), result.end(), [](int p) { return p >= 3 && p <= 10; }));
}

TEST(StochasticRewriteTest, TestCalculateRs) {
  const std::vector<size_t> partition{3, 4, 3, 5};

  const std::vector<size_t> rs{1, 4, 8, 11, 16};

  EXPECT_EQ(rs, calculateRs(partition));
}


TEST(StochasticRewriteTest, TestCalculateYGens) {
  const std::vector<size_t> partition{3, 4, 3, 5};

  const auto y_gens = calculateYGensInBGens(partition);

  ASSERT_EQ(15, y_gens.size());

  EXPECT_EQ(y_gens[0], std::vector<int>({1, 2, 3}));
  EXPECT_EQ(y_gens[1], std::vector<int>({2, 3}));
  EXPECT_EQ(y_gens[2], std::vector<int>({3}));
  EXPECT_EQ(y_gens[3], std::vector<int>({4, 5, 6, 7}));
  EXPECT_EQ(y_gens[4], std::vector<int>({5, 6, 7}));
  EXPECT_EQ(y_gens[5], std::vector<int>({6, 7}));
  EXPECT_EQ(y_gens[6], std::vector<int>({7}));
  EXPECT_EQ(y_gens[7], std::vector<int>({8, 9, 10}));
  EXPECT_EQ(y_gens[8], std::vector<int>({9, 10}));
  EXPECT_EQ(y_gens[9], std::vector<int>({10}));
  EXPECT_EQ(y_gens[10], std::vector<int>({11, 12, 13, 14, 15}));
  EXPECT_EQ(y_gens[11], std::vector<int>({12, 13, 14, 15}));
  EXPECT_EQ(y_gens[12], std::vector<int>({13, 14, 15}));
  EXPECT_EQ(y_gens[13], std::vector<int>({14, 15}));
  EXPECT_EQ(y_gens[14], std::vector<int>({15}));
}


TEST(StochasticRewriteTest, TestCalculateBGens) {
  const std::vector<size_t> partition{3, 4, 3, 5};

  const auto b_gens = calculateBGensInYGens(partition);

  ASSERT_EQ(15, b_gens.size());

  EXPECT_EQ(b_gens[0], std::vector<int>({1, -2}));
  EXPECT_EQ(b_gens[1], std::vector<int>({2, -3}));
  EXPECT_EQ(b_gens[2], std::vector<int>({3}));
  EXPECT_EQ(b_gens[3], std::vector<int>({4, -5}));
  EXPECT_EQ(b_gens[4], std::vector<int>({5, -6}));
  EXPECT_EQ(b_gens[5], std::vector<int>({6, -7}));
  EXPECT_EQ(b_gens[6], std::vector<int>({7}));
  EXPECT_EQ(b_gens[7], std::vector<int>({8, -9}));
  EXPECT_EQ(b_gens[8], std::vector<int>({9, -10}));
  EXPECT_EQ(b_gens[9], std::vector<int>({10}));
  EXPECT_EQ(b_gens[10], std::vector<int>({11, -12}));
  EXPECT_EQ(b_gens[11], std::vector<int>({12, -13}));
  EXPECT_EQ(b_gens[12], std::vector<int>({13, -14}));
  EXPECT_EQ(b_gens[13], std::vector<int>({14, -15}));
  EXPECT_EQ(b_gens[14], std::vector<int>({15}));
}


TEST(StochasticRewriteTest, TestYRelations) {
  const std::vector<size_t> partition{3};

  const auto y_rels = calculateYRelations(partition);

  ASSERT_EQ(3, y_rels.size());

  EXPECT_EQ(y_rels[0], std::vector<int>({1, -3, 1, -2, 3, -1, 3, -2}));
  EXPECT_EQ(y_rels[1], std::vector<int>({1, -2, 3, 2, -1, -3}));
  EXPECT_EQ(y_rels[2], std::vector<int>({2, 2, -3, -2, -3}));
}


TEST(StochasticRewriteTest, TestAdditionalRelations) {
  const std::vector<size_t> partition{4, 4};

  const auto a_rels = calculateAdditionalRelations(partition);

  ASSERT_EQ(12, a_rels.size());
  EXPECT_EQ(a_rels[0], std::vector<int>({2, 1, 4, -1, -1}));
  EXPECT_EQ(a_rels[1], std::vector<int>({3, 1, 4, -2, -1}));
  EXPECT_EQ(a_rels[2], std::vector<int>({4, 1, 4, -3, -1}));
  EXPECT_EQ(a_rels[3], std::vector<int>({3, 2, 4, -2, -2}));
  EXPECT_EQ(a_rels[4], std::vector<int>({4, 2, 4, -3, -2}));
  EXPECT_EQ(a_rels[5], std::vector<int>({4, 3, 4, -3, -3}));
  EXPECT_EQ(a_rels[6], std::vector<int>({6, 5, 8, -5, -5}));
  EXPECT_EQ(a_rels[7], std::vector<int>({7, 5, 8, -6, -5}));
  EXPECT_EQ(a_rels[8], std::vector<int>({8, 5, 8, -7, -5}));
  EXPECT_EQ(a_rels[9], std::vector<int>({7, 6, 8, -6, -6}));
  EXPECT_EQ(a_rels[10], std::vector<int>({8, 6, 8, -7, -6}));
  EXPECT_EQ(a_rels[11], std::vector<int>({8, 7, 8, -7, -7}));
}


TEST(StochasticRewriteTest, TestSRelations) {
  const std::vector<size_t> partition{3};

  const auto s_rels = calculateSRelations(partition);

  ASSERT_EQ(34, s_rels.size());
}


TEST(StochasticRewriteTest, TestRules) {
  const std::vector<size_t> partition{3, 4};

  const auto rules = calculateRewritingRules(partition);

  ASSERT_EQ(380, rules.size());
}


TEST(StochasticRewriteTest, TestSelectSubwords) {
  const std::vector<size_t> partition{10, 13, 4, 3, 5};

  std::mt19937 g(1234);

  const auto rs = calculateRs(partition, 0);

  for (size_t i = 0; i < 10; ++i) {
    const auto result = selectSubwords(partition, g);

    ASSERT_EQ(5, result.size());

    for (size_t i = 0; i < result.size(); ++i) {
      EXPECT_TRUE(result[i] >= rs[i] && result[i] + 2 <= rs[i + 1]);
    }
  }
}


TEST(StochasticRewriteTest, TestRewrite) {
  const std::vector<size_t> partition{3, 4};

  const std::vector<int> w{1, 2, 1, 2, 3};

  std::mt19937 g(1234);

  const auto rules = calculateRewritingRules(partition);

  const auto rewritten = rewrite(w, 3, 6, rules, g);

  EXPECT_EQ(std::vector<int>({3, 1, 3, 1, 2, 3}), rewritten);
}


static bool compareBraids(int n, const Word& lhs, const Word& rhs) {
  LinkedBraidStructure lbs(n - 1, lhs * -rhs);

  lbs.removeLeftHandles();

  return lbs.size() == 0;
}


TEST(StochasticRewriteTest, TestStochasticRewrite1) {
  const std::vector<size_t> partition{3, 4};

  const Word w({1, 2, 1, 2, 3});

  std::mt19937 g(1234);

  const auto rewritten = stochasticRewrite(w, partition, 3, 6, 3, g);

  EXPECT_TRUE(compareBraids(7, w, rewritten));
}


TEST(StochasticRewriteTest, TestStochasticRewrite2) {
  const size_t n = 16;

  std::mt19937 g(1234);

  for (size_t i = 0; i < 10; ++i) {
    const auto partition = randomPartition(n - 1, 3, n - 1, g);

    const auto w = random::randomWord(n - 1, 1000, g);

    const auto rewritten = stochasticRewrite(w, partition, 5, 10, 3, g);

    EXPECT_TRUE(compareBraids(n, w, rewritten));
  }
}


TEST(StochasticRewriteTest, TestStochasticRewrite3) {
  const size_t n = 16;

  std::mt19937 g(1234);

  const auto partition = randomPartition(n - 1, 3, n - 1, g);

  const auto w = random::randomWord(n - 1, 20000, g);

  const auto rewritten = stochasticRewrite(w, partition, 5, 10, 3, g);

  EXPECT_TRUE(compareBraids(n, w, rewritten));
}
} // namespace
} // namespace stochasticrewrite
} // namespace crag
