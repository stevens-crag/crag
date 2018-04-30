#include "gtest/gtest.h"

#include "WordRep.h"

namespace {

TEST(WordRep, RequireNonzeroIndices) {
  EXPECT_THROW({ WordRep(0); }, std::invalid_argument);
  EXPECT_THROW({ WordRep({1, 2, -1, 0}); }, std::invalid_argument);
}

TEST(WordRep, Print) {
  EXPECT_EQ("x1 x2 x3", WordRep({1, 2, 3}).toString());
  EXPECT_EQ("x2", WordRep({1, -1, 2}).toString());
  EXPECT_EQ("1", WordRep({1, -1}).toString());
  EXPECT_EQ("1", WordRep().toString());
  EXPECT_EQ("1", WordRep({1, 1, -1, -1}).toString());
}

TEST(WordRep, Basics) {
  const WordRep w({1, 2, 3});

  ASSERT_EQ(3, w.size());
  ASSERT_EQ(3, w.length());
  ASSERT_FALSE(w.empty());

  auto it = w.begin();

  EXPECT_EQ(1, *(it++));
  EXPECT_EQ(2, *(it++));
  EXPECT_EQ(3, *(it++));
  EXPECT_EQ(w.end(), it);

  auto rit = w.rbegin();

  EXPECT_EQ(3, *(rit++));
  EXPECT_EQ(2, *(rit++));
  EXPECT_EQ(1, *(rit++));
  EXPECT_EQ(w.rend(), rit);

  auto u = w;

  u.clear();

  EXPECT_TRUE(u.empty());
}

TEST(WordRep, RelationalOperators) {
  const WordRep u({1, 2, 3});
  const WordRep v({1, 2, 3, 4});

  EXPECT_EQ(WordRep({1, 2, 3}), u);
  EXPECT_NE(u, v);

  EXPECT_LT(u, v);
}

TEST(WordRep, Conjugation) {
  EXPECT_EQ(WordRep({2, 1, 1}), WordRep({1, 2, 1})^WordRep({1}));
  EXPECT_EQ(WordRep({-2, 1}), WordRep({1, -2})^WordRep({2}));
  EXPECT_EQ(WordRep({-3, 1, 2, 3}), WordRep({1, 2})^WordRep({3}));
}

TEST(WordRep, Power) {
  EXPECT_EQ(WordRep(), WordRep({2, 1, 3})^0);
  EXPECT_EQ(WordRep({2, 1, 2, 1}), WordRep({2, 1})^2);
  EXPECT_EQ(WordRep({-1, 2, 2, 2, 1}), WordRep({-1, 2, 1})^3);
  EXPECT_EQ(WordRep({-1, -2, -2, -2, 1}), WordRep({-1, 2, 1})^-3);
}

TEST(WordRep, Multiplication) {
  EXPECT_EQ(WordRep({1, 2}), WordRep({1, 2, -1}) * WordRep({1}));
  EXPECT_EQ(WordRep({1, 2, 1, 3}), WordRep({1, 2}) * WordRep({1, 3}));
}

TEST(WordRep, Inverse) {
  EXPECT_EQ(WordRep({1, -2, -1}), WordRep({1, 2, -1}).inverse());
}

TEST(WordRep, ExponentSums) {
  const WordRep w({1, 1, 2, 3, -1, -2});

  EXPECT_EQ(1, w.exponentSum(1));
  EXPECT_EQ(0, w.exponentSum(2));
  EXPECT_EQ(1, w.exponentSum(3));

  EXPECT_THROW({ w.exponentSum(0); }, std::invalid_argument);
}

TEST(WordRep, Occurrences) {
  const WordRep w({1, 1, 2, 3, -1, -2});

  EXPECT_EQ(3, w.occurrences(1));
  EXPECT_EQ(2, w.occurrences(2));
  EXPECT_EQ(1, w.occurrences(3));

  EXPECT_THROW({ w.occurrences(0); }, std::invalid_argument);
}

TEST(WordRep, Contains) {
  EXPECT_TRUE(WordRep(1).contains(1));
  EXPECT_TRUE(WordRep(-1).contains(1));
  EXPECT_FALSE(WordRep({1, 3}).contains(2));
}

TEST(WordRep, Root_01) {
  WordRep w({1, 2, 1, 2});
  const auto p = w.root();

  EXPECT_EQ(WordRep({1, 2}), p.first);
  EXPECT_EQ(2, p.second);
}

TEST(WordRep, DISABLED_Root_02) {
  WordRep w({-3, -4, 1, 2, 1, 2, 4, 3});
  const auto p = w.root();

  EXPECT_EQ(WordRep({-3, -4, 1, 2, 4, 3}), p.first);
  EXPECT_EQ(2, p.second);
}

TEST(WordRep, CyclicallyReduce_01) {
  WordRep w({-1, 2, 3, -2, 1});
  const auto conjugator = w.cyclicallyReduce();

  EXPECT_EQ(WordRep(3), w);
  EXPECT_EQ(WordRep({-2, 1}), conjugator);
}

TEST(WordRep, CyclicallyReduce_02) {
  WordRep w({1, 2, 3});
  const auto conjugator = w.cyclicallyReduce();

  EXPECT_EQ(WordRep({1, 2, 3}), w);
  EXPECT_EQ(WordRep(), conjugator);
}

TEST(WordRep, CyclicLeftShift_01) {
  WordRep w({1, 2, 3});

  w.cyclicLeftShift();
  EXPECT_EQ(WordRep({2, 3, 1}), w);

  w.cyclicLeftShift();
  EXPECT_EQ(WordRep({3, 1, 2}), w);

  w.cyclicLeftShift();
  EXPECT_EQ(WordRep({1, 2, 3}), w);
}

TEST(WordRep, CyclicLeftShift_02) {
  WordRep w({1, 2, 3, -2, -1});

  w.cyclicLeftShift();
  EXPECT_EQ(WordRep({2, 3, -2}), w);

  w.cyclicLeftShift();
  EXPECT_EQ(WordRep({3}), w);
}

TEST(WordRep, CyclicRightShift_01) {
  WordRep w({1, 2, 3});

  w.cyclicRightShift();
  EXPECT_EQ(WordRep({3, 1, 2}), w);

  w.cyclicRightShift();
  EXPECT_EQ(WordRep({2, 3, 1}), w);

  w.cyclicRightShift();
  EXPECT_EQ(WordRep({1, 2, 3}), w);
}

TEST(WordRep, CyclicRightShift_02) {
  WordRep w({1, 2, 3, -2, -1});

  w.cyclicRightShift();
  EXPECT_EQ(WordRep({2, 3, -2}), w);

  w.cyclicRightShift();
  EXPECT_EQ(WordRep({3}), w);
}

TEST(WordRep, CyclicPermutation_01) {
  WordRep w({1, 2, 3, -2, -1});
  w.cyclicallyPermute(2);
  EXPECT_EQ(WordRep(3), w);

  w = WordRep({1, 2, 3, -2, -1});
  w.cyclicallyPermute(4);
  EXPECT_EQ(WordRep({2, 3, -2}), w);
}

TEST(WordRep, CyclicPermutation_02) {
  WordRep w({1, 2, 3, -2, -1});
  w.cyclicallyPermute(-2);
  EXPECT_EQ(WordRep(3), w);
}

TEST(WordRep, InitialSegment) {
  WordRep w({1, 2, 3, 4, 5, 6, 7});
  w.initialSegment(3);
  EXPECT_EQ(WordRep({1, 2, 3}), w);

  w.initialSegment(0);
  EXPECT_EQ(WordRep(), w);

  w = WordRep({1, 2, 3, 4, 5, 6, 7});
  w.initialSegment(10);
  EXPECT_EQ(WordRep({1, 2, 3, 4, 5, 6, 7}), w);

  w = WordRep({1, 2, 3, 4, 5, 6, 7});
  w.initialSegment(7);
  EXPECT_EQ(WordRep({1, 2, 3, 4, 5, 6, 7}), w);
}

TEST(WordRep, TerminalSegment) {
  WordRep w({1, 2, 3, 4, 5, 6, 7});
  w.terminalSegment(5);
  EXPECT_EQ(WordRep({6, 7}), w);

  w.terminalSegment(100);
  EXPECT_EQ(WordRep(), w);

  w = WordRep({1, 2, 3, 4, 5, 6, 7});
  w.terminalSegment(0);
  EXPECT_EQ(WordRep({1, 2, 3, 4, 5, 6, 7}), w);
}

TEST(WordRep, Segment) {
  WordRep w({1, 2, 3, 4, 5, 6, 7});
  EXPECT_EQ(WordRep(), w.subword(0, 0));
  w.segment(0, 0);
  EXPECT_EQ(WordRep(), w);

  w = WordRep({1, 2, 3, 4, 5, 6, 7});
  EXPECT_EQ(WordRep(), w.subword(7, 10));
  w.segment(7, 10);
  EXPECT_EQ(WordRep(), w);

  w = WordRep({1, 2, 3, 4, 5, 6, 7});
  EXPECT_EQ(WordRep({6, 7}), w.subword(5, 100));
  w.segment(5, 100);
  EXPECT_EQ(WordRep({6, 7}), w);

  w = WordRep({1, 2, 3, 4, 5, 6, 7});
  EXPECT_EQ(WordRep({2, 3}), w.subword(1, 3));
  w.segment(1, 3);
  EXPECT_EQ(WordRep({2, 3}), w);
}

TEST(WordRep, Insert_01) {
  WordRep w({1, 2});

  w.insert(1, 3);
  EXPECT_EQ(WordRep({1, 3, 2}), w);

  w.insert(w.end(), 4);
  EXPECT_EQ(WordRep({1, 3, 2, 4}), w);

  w.insert(100, 5);
  EXPECT_EQ(WordRep({1, 3, 2, 4, 5}), w);

  w.insert(w.begin(), 3);
  EXPECT_EQ(WordRep({3, 1, 3, 2, 4, 5}), w);

  EXPECT_THROW({ w.insert(w.begin(), 0); }, std::invalid_argument);
}

TEST(WordRep, Insert_02) {
  WordRep w({1, 2, 3});

  w.insert(w.end(), -3);
  EXPECT_EQ(std::list<int>({1, 2, 3, -3}), w.toList());

  w.freelyReduce();
  EXPECT_EQ(WordRep({1, 2}), w);

  w.insert(1, -1);
  EXPECT_EQ(std::list<int>({1, -1, 2}), w.toList());

  w.freelyReduce();
  EXPECT_EQ(WordRep(2), w);
}

TEST(WordRep, Insert_03) {
  WordRep w({1, 2, 3});

  w.insert(w.begin(), -1);
  EXPECT_EQ(std::list<int>({-1, 1, 2, 3}), w.toList());

  w.freelyReduce();
  EXPECT_EQ(WordRep({2, 3}), w);

  w.insert(1, -3);
  EXPECT_EQ(std::list<int>({2, -3, 3}), w.toList());

  w.freelyReduce();
  EXPECT_EQ(WordRep(2), w);
}

TEST(WordRep, InsertRange_01) {
  WordRep w({1, 2, 3});

  std::vector<int> v = {1, 2};

  w.insert(w.end(), v.begin(), v.end());

  EXPECT_EQ(WordRep({1, 2, 3, 1, 2}), w);

  std::vector<int> zero_vec(5, 0);

  EXPECT_THROW({ w.insert(w.begin(), zero_vec.begin(), zero_vec.end()); }, std::invalid_argument);
}

TEST(WordRep, InsertRange_02) {
  WordRep w({1, 2});

  std::vector<int> v = {1, 2};

  w.insert(w.begin(), v.begin(), v.end());
  EXPECT_EQ(WordRep({1, 2, 1, 2}), w);
}

TEST(WordRep, InsertRange_03) {
  WordRep w({4, 1, 2, 3});

  std::vector<int> v = {-1, -2};

  w.insert(2, v.begin(), v.end());
  EXPECT_EQ(std::list<int>({4, 1, -1, -2, 2, 3}), w.toList());

  w.freelyReduce(++w.begin(), --w.end());
  EXPECT_EQ(WordRep({4, 3}), w);
}

TEST(WordRep, Replace) {
  WordRep w({1, 2, 3});
  w.replace(1, 3);
  EXPECT_EQ(WordRep({1, 3, 3}), w);

  w.replace(++w.begin(), 4);
  EXPECT_EQ(WordRep({1, 4, 3}), w);

  w.replace(0, -4);
  EXPECT_EQ(std::list<int>({-4, 4, 3}), w.toList());

  w.freelyReduce(w.begin(), --w.end());
  EXPECT_EQ(WordRep(3), w);

  EXPECT_THROW({ w.replace(w.begin(), 0); }, std::invalid_argument);
  EXPECT_THROW({ w.replace(100, 1); }, std::invalid_argument);
}

TEST(WordRep, ReplaceRange) {
  WordRep w({1, 2, 3});

  std::vector<int> v = {6, 7, 8};

  w.replace(++w.begin(), v.begin(), v.end());

  EXPECT_EQ(WordRep({1, 6, 7}), w);

  w.replace(0, v.begin(), v.end());

  EXPECT_EQ(WordRep({6, 7, 8}), w);

  EXPECT_THROW({ w.replace(100, v.begin(), v.end()); }, std::invalid_argument);
}
} // namespace
