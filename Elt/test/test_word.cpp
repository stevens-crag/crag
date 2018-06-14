#include "gtest/gtest.h"

#include "Word.h"

namespace {

TEST(Word, TestInput1) {
  const auto w = "x1 x2^-2"_w;
  const Word result({1, -2, -2});

  EXPECT_EQ(result, w);
}

TEST(Word, Print) {
  EXPECT_EQ("x1 x2 x3", Word({1, 2, 3}).toString());
  EXPECT_EQ("x2", Word({1, -1, 2}).toString());
  EXPECT_EQ("1", Word({1, -1}).toString());
  EXPECT_EQ("1", Word().toString());
  EXPECT_EQ("1", Word({1, 1, -1, -1}).toString());
}

TEST(Word, Basics) {
  const Word w({1, 2, 3});

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

  EXPECT_EQ(std::list<int>({1, 2, 3}), w.toList());
}

TEST(Word, Conjugation_01) {
  Word w({1, 2});
  w ^= Word({3});
  EXPECT_EQ(Word({-3, 1, 2, 3}), w);
}

TEST(Word, Conjugation_02) {
  Word w({1, 2});
  auto u = w ^ Word({3});

  EXPECT_EQ(Word({1, 2}), w);
  EXPECT_EQ(Word({-3, 1, 2, 3}), u);
}

TEST(Word, Power_01) {
  Word w({1, 2});
  w ^= 2;
  EXPECT_EQ(Word({1, 2, 1, 2}), w);
}

TEST(Word, Power_02) {
  Word w({1, 2});
  auto u = w ^ 2;

  EXPECT_EQ(Word({1, 2}), w);
  EXPECT_EQ(Word({1, 2, 1, 2}), u);
}

TEST(Word, Power_03) {
  Word w({1, 2});

  auto u = w.power(0);
  EXPECT_EQ(Word({1, 2}), w);
  EXPECT_EQ(Word(), u);

  u = w.power(2);
  EXPECT_EQ(Word({1, 2}), w);
  EXPECT_EQ(Word({1, 2, 1, 2}), u);

  u = w.power(-1);
  EXPECT_EQ(Word({1, 2}), w);
  EXPECT_EQ(Word({-2, -1}), u);
}

TEST(Word, Multiplication_01) {
  Word w({1, 2});
  w *= Word({3, 4});
  EXPECT_EQ(Word({1, 2, 3, 4}), w);
}

TEST(Word, Multiplication_02) {
  Word w({1, 2});
  auto u = w * Word({3, 4});

  EXPECT_EQ(Word({1, 2}), w);
  EXPECT_EQ(Word({1, 2, 3, 4}), u);
}

TEST(Word, Inverse) {
  Word w({1, 2});

  EXPECT_EQ(Word({-2, -1}), -w);
  EXPECT_EQ(Word({-2, -1}), w.inverse());
}

TEST(Word, PushBack_01) {
  Word w({1, 2});

  const auto u = w;

  w.push_back(-2);

  EXPECT_EQ(Word(1), w);
  EXPECT_EQ(Word({1, 2}), u);
}

TEST(Word, PushBack_02) {
  Word w({1, 2});

  const auto u = w;

  w.push_back(Word({-2, 1, 3}));

  EXPECT_EQ(Word({1, 1, 3}), w);
  EXPECT_EQ(Word({1, 2}), u);
}

TEST(Word, PushFront_01) {
  Word w({1, 2});

  const auto u = w;

  w.push_front(-1);

  EXPECT_EQ(Word(2), w);
  EXPECT_EQ(Word({1, 2}), u);
}

TEST(Word, PushFront_02) {
  Word w({1, 2});

  const auto u = w;

  w.push_front(Word({3, 2, -1}));

  EXPECT_EQ(Word({3, 2, 2}), w);
  EXPECT_EQ(Word({1, 2}), u);
}

TEST(Word, Pop) {
  Word w({1, 2, 3});

  w.pop_back();
  EXPECT_EQ(Word({1, 2}), w);

  w.pop_front();
  EXPECT_EQ(Word(2), w);
}

TEST(Word, Root) {
  Word w({1, 2, 1, 2, 1, 2});
  Word base;

  const auto exp = w.getPower(base);

  EXPECT_EQ(Word({1, 2}), base);
  EXPECT_EQ(3, exp);
}

TEST(Word, DoesContain) {
  Word w({1, 2});

  EXPECT_TRUE(w.doesContain(1));
  EXPECT_TRUE(w.doesContain(-2));

  EXPECT_FALSE(w.doesContain(3));

  EXPECT_THROW({ w.doesContain(0); }, std::invalid_argument);
}

TEST(Word, CyclicallyReduce_01) {
  Word w({-3, 1, 2, 3});

  const auto w_reduced = w.cyclicallyReduce();
  EXPECT_EQ(Word({1, 2}), w_reduced);

  w.cyclicallyReduceWord();
  EXPECT_EQ(w_reduced, w);
}

TEST(Word, CyclicallyReduce_02) {
  Word w({-3, 1, 2, 3});

  Word c;

  const auto w_reduced = w.cyclicallyReduce(c);
  EXPECT_EQ(Word({1, 2}), w_reduced);
  EXPECT_EQ(Word(3), c);

  c = Word();
  EXPECT_EQ(Word(), c);

  w.cyclicallyReduceWord(c);
  EXPECT_EQ(w_reduced, w);
  EXPECT_EQ(Word(3), c);
}

TEST(Word, ExponentSums) {
  const Word w({1, 1, 2, 3, -1, -2});

  EXPECT_EQ(1, w.exponentSum(1));
  EXPECT_EQ(0, w.exponentSum(2));
  EXPECT_EQ(1, w.exponentSum(3));

  EXPECT_THROW({ w.exponentSum(0); }, std::invalid_argument);
}

TEST(Word, Occurrences_1) {
  const Word w({1, 1, 2, 3, -1, -2});

  EXPECT_EQ(3, w.isIn(1));
  EXPECT_EQ(2, w.isIn(2));
  EXPECT_EQ(1, w.isIn(3));

  EXPECT_THROW({ w.isIn(0); }, std::invalid_argument);
}

TEST(Word, Occurrences_2) {
  using map = std::map<size_t, size_t>;

  EXPECT_EQ(map({
      {1, 1},
  }), occurrences(Word({1})));

  EXPECT_EQ(map({
      {1, 2},
      {2, 1},
  }), occurrences(Word({1, 2, -1})));

  EXPECT_EQ(map({
      {1, 2},
      {2, 2},
      {3, 1},
      {4, 2},
  }), occurrences(Word({1, 2, 4, -2, -4, -1, 3})));
}

TEST(Word, Insert_01) {
  Word w({1, 2});
  auto u = w;

  w.insert(1, -2);
  EXPECT_EQ(std::list<int>({1, -2, 2}), w.toList());
  EXPECT_EQ(Word({1, 2}), u);

  w = w.freelyReduce();
  EXPECT_EQ(Word(1), w);
}

TEST(Word, Insert_02) {
  Word w({1, 2});
  auto u = w;

  w.insert(2, 3);
  EXPECT_EQ(Word({1, 2, 3}), w);
  EXPECT_EQ(Word({1, 2}), u);

  w.insert(100, 4);
  EXPECT_EQ(Word({1, 2, 3, 4}), w);

  EXPECT_THROW({ w.insert(0, 0); }, std::invalid_argument);
}

TEST(Word, InsertRange_01) {
  Word w({1, 2, 3});

  std::vector<int> v = {1, 2};

  w.insert(3, v.begin(), v.end());

  EXPECT_EQ(Word({1, 2, 3, 1, 2}), w);

  std::vector<int> zero_vec(5, 0);

  EXPECT_THROW({ w.insert(0, zero_vec.begin(), zero_vec.end()); }, std::invalid_argument);
}

TEST(Word, InsertRange_02) {
  Word w({4, 1, 2, 3});

  std::vector<int> v = {-1, -2};

  w.insert(2, v.begin(), v.end());
  EXPECT_EQ(std::list<int>({4, 1, -1, -2, 2, 3}), w.toList());

  w = w.freelyReduce(++w.begin(), --w.end());
  EXPECT_EQ(Word({4, 3}), w);
}

TEST(Word, Replace) {
  Word w({1, 2, 3});
  w.replace(1, 3);
  EXPECT_EQ(Word({1, 3, 3}), w);

  w.replace(1, 4);
  EXPECT_EQ(Word({1, 4, 3}), w);

  w.replace(0, -4);
  EXPECT_EQ(std::list<int>({-4, 4, 3}), w.toList());

  w = w.freelyReduce(w.begin(), --w.end());
  EXPECT_EQ(Word(3), w);

  EXPECT_THROW({ w.replace(0, 0); }, std::invalid_argument);
  EXPECT_THROW({ w.replace(100, 1); }, std::invalid_argument);
}

TEST(Word, ReplaceRange) {
  Word w({1, 2, 3});

  std::vector<int> v = {6, 7, 8};

  w.replace(1, v.begin(), v.end());

  EXPECT_EQ(Word({1, 6, 7}), w);

  w.replace(0, v.begin(), v.end());

  EXPECT_EQ(Word({6, 7, 8}), w);

  EXPECT_THROW({ w.replace(100, v.begin(), v.end()); }, std::invalid_argument);
}

TEST(Word, ReplaceGenerators_01) {
  Word w({1, 2});

  std::vector<Word> map{
      Word({2, -1}),
      Word({1, 2, 1, 1}),
  };

  const auto image = w.replaceGenerators(map);

  EXPECT_EQ(Word({2, 2, 1, 1}), image);
}

TEST(Word, ReplaceGenerators_02) {
  Word w({1, 3, 2});

  std::vector<Word> map{
      Word({2, -1}),
      Word({1, 2, 1, 1}),
  };

  EXPECT_THROW({ w.replaceGenerators(map); }, std::invalid_argument);
}

TEST(Word, Abelianization) {
  EXPECT_EQ(Word(), abelianization(Word()));
  EXPECT_EQ(Word(1), abelianization(Word(1)));
  EXPECT_EQ(Word({1, 2}), abelianization(Word({2, 1})));
  EXPECT_EQ(Word({1, 3, 3}), abelianization(Word({2, 3, -2, 3, 1})));
  EXPECT_EQ(Word({1, 1, 1, 2, 2, 2}), abelianization(Word({1, 2, 1, 2, 1, 2})));
  EXPECT_EQ(Word(), abelianization(Word({1, 2, -1, -2})));
}
} // namespace
