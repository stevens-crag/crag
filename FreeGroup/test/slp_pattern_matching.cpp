/*
 * SLPSet_inspector.cpp
 *
 *  Created on: Jan 22, 2013
 *      Author: dpantele
 */

#include <memory>
#include <functional>
#include <utility>
#include <vector>
#include <tuple>
#include <algorithm>

#include "gtest/gtest.h"
#include "slp_inspector.h"
#include "slp_pattern_matching.h"
#include "slp_vertex_word.h"

namespace crag {
namespace slp {
namespace {
typedef ::std::tuple<int, int, int, int, ::std::vector<int>> BoundedTaskAcceptorTestParam;

class BoundedTaskAcceptorTest : public ::testing::TestWithParam<BoundedTaskAcceptorTestParam>
{  };

TEST_P(BoundedTaskAcceptorTest, CheckOrder) {
  int pattern_code;
  int text_code;
  int lookup_from_int;
  int lookup_length_int;
  ::std::vector<int> correct_path;

  ::std::tie(pattern_code, text_code, lookup_from_int, lookup_length_int, correct_path) = GetParam();

  TerminalVertex t(TerminalSymbol{} + 1);
  NonterminalVertex t2(t, t);
  NonterminalVertex t4(t2, t2);

  Vertex pattern = Vertex(), text = Vertex();
  switch(pattern_code) {
  case 1:
    pattern = t;
    break;
  case 2:
    pattern = t2;
    break;
  }

  switch(text_code) {
  case 1:
    text = t;
    break;
  case 2:
    text = t2;
    break;
  case 4:
    text = t4;
    break;
  }

  //Must update this code if you change constructors in PatternMatchesGenerator
  LongInteger lookup_from = lookup_from_int;
  LongInteger lookup_length = lookup_length_int;
  LongInteger first_lookup_begin_position_ = lookup_from;
  LongInteger first_lookup_end_position_(lookup_from + pattern.length());
  LongInteger last_lookup_begin_position_(((lookup_from += lookup_length) > text.length()) ? (text.length() - pattern.length()) : (lookup_from - pattern.length()));
  Inspector<inspector::Inorder, inspector::BoundedTaskAcceptor> text_inspector_(text, inspector::BoundedTaskAcceptor(
    first_lookup_end_position_, last_lookup_begin_position_, pattern.length()
  ));

  Inspector<inspector::Inorder, inspector::BoundedTaskAcceptor>& text_inspector = text_inspector_;

  for (auto correct_vertex : correct_path) {
    ASSERT_FALSE(text_inspector.stopped());
    if (correct_vertex == 1) {
      EXPECT_EQ(t, text_inspector.vertex());
    } else if (correct_vertex == 2) {
      EXPECT_EQ(t2, text_inspector.vertex());
    } else if (correct_vertex == 4) {
      EXPECT_EQ(t4, text_inspector.vertex());
    } else {
      FAIL() << "Wrong vertex id in vertices order vector";
    }

    text_inspector.next();
  }
  ASSERT_TRUE(text_inspector.stopped());
}

INSTANTIATE_TEST_CASE_P(
    BoundedTaskAcceptor,
    BoundedTaskAcceptorTest,
    ::testing::Values(
      ::std::make_tuple(1, 0, 0, 1, ::std::vector<int>()),
      ::std::make_tuple(1, 1, 0, 1, ::std::vector<int>({1})),
      ::std::make_tuple(2, 1, 0, 2, ::std::vector<int>({})),
      ::std::make_tuple(2, 1, 0, 1, ::std::vector<int>({})),
      ::std::make_tuple(1, 1, 1, 1, ::std::vector<int>({})),
      ::std::make_tuple(1, 2, 0, 10, ::std::vector<int>({1, 2, 1})),
      ::std::make_tuple(2, 2, 0, 10, ::std::vector<int>({2})),
      ::std::make_tuple(1, 2, 0, 1, ::std::vector<int>({1, 2})),
      ::std::make_tuple(1, 2, 1, 1, ::std::vector<int>({2, 1})),
      ::std::make_tuple(2, 4, 0, 4, ::std::vector<int>({2, 4, 2}))
    ));

class FilledMatchingTable : public MatchingTable {
  public:
    FilledMatchingTable(const Vertex& pattern, const Vertex& text)
      : MatchingTable()
    {
      PostorderInspector text_inspector(text);

      while (!text_inspector.stopped()) {
        VertexWord current_text_word(text_inspector.vertex());
        std::string current_text(current_text_word.begin(), current_text_word.end());

        PostorderInspector pattern_inspector(pattern);
        while (!pattern_inspector.stopped()) {
          FiniteArithmeticSequence result;

          if (pattern_inspector.vertex().length() <= text_inspector.vertex().length() &&
              match_table_->find(::std::make_pair(pattern_inspector.vertex(), text_inspector.vertex())) == match_table_->end()) {
            VertexWord current_pattern_word(pattern_inspector.vertex());
            std::string current_pattern(current_pattern_word.begin(), current_pattern_word.end());
            size_t current_match = text_inspector.vertex().split_point().get_ui() -
                pattern_inspector.vertex().length().get_ui();

            if (text_inspector.vertex().split_point() < pattern_inspector.vertex().length()) {
              current_match = 0;
            }

            size_t last_possible_match = text_inspector.vertex().split_point().get_ui();
//            std::cout << get_vertex_text_as_string(pattern_inspector.current_vertex()) << " in "
//                << get_vertex_text_as_string(text_inspector.current_vertex())
//                << " from " << current_match << " to " << last_possible_match << std::endl;

            size_t first_match = current_text.find(current_pattern, current_match);
            size_t count = 0;
            current_match = first_match;
            size_t step = 0;
            while (current_match <= last_possible_match && current_match != std::string::npos) {
              if (step == 0) {
                step = current_match - first_match;
              }
              ++count;
              ++current_match;
              current_match = current_text.find(current_pattern, current_match);
            }
            if (first_match != ::std::string::npos) {
              if (step == 0) {
                step = 1;
              }
              result = FiniteArithmeticSequence(first_match, step, count);
            }
          }
          this->match_table_->insert(::std::make_pair(::std::make_pair(pattern_inspector.vertex(), text_inspector.vertex()), result));
          pattern_inspector.next();
        }
        text_inspector.next();
      }
    }
};

TEST(FilledMatchingTable, Example1) {
  TerminalVertex a('a');
  NonterminalVertex a2(a, a);
  NonterminalVertex a4(a2, a2);
  NonterminalVertex text(a4, a4);

  NonterminalVertex a3(a, a2);
  NonterminalVertex pattern(a3, a2);

  FilledMatchingTable matching_table(pattern, text);

  EXPECT_EQ(FiniteArithmeticSequence(0, 1, 1), matching_table.matches(a,a));
  EXPECT_EQ(FiniteArithmeticSequence(), matching_table.matches(a2,a));
  EXPECT_EQ(FiniteArithmeticSequence(), matching_table.matches(a3,a));
  EXPECT_EQ(FiniteArithmeticSequence(), matching_table.matches(pattern,a));
  EXPECT_EQ(FiniteArithmeticSequence(0, 1, 2), matching_table.matches(a,a2));
  EXPECT_EQ(FiniteArithmeticSequence(0, 1, 1), matching_table.matches(a2,a2));
  EXPECT_EQ(FiniteArithmeticSequence(), matching_table.matches(a3,a2));
  EXPECT_EQ(FiniteArithmeticSequence(), matching_table.matches(pattern,a2));
  EXPECT_EQ(FiniteArithmeticSequence(1, 1, 2), matching_table.matches(a,a4));
  EXPECT_EQ(FiniteArithmeticSequence(0, 1, 3), matching_table.matches(a2,a4));
  EXPECT_EQ(FiniteArithmeticSequence(0, 1, 2), matching_table.matches(a3,a4));
  EXPECT_EQ(FiniteArithmeticSequence(), matching_table.matches(pattern,a4));
  EXPECT_EQ(FiniteArithmeticSequence(3, 1, 2), matching_table.matches(a,text));
  EXPECT_EQ(FiniteArithmeticSequence(2, 1, 3), matching_table.matches(a2,text));
  EXPECT_EQ(FiniteArithmeticSequence(1, 1, 4), matching_table.matches(a3,text));
  EXPECT_EQ(FiniteArithmeticSequence(0, 1, 4), matching_table.matches(pattern,text));
}

TEST(FilledMatchingTable, Example2) {
  TerminalVertex a('a');
  TerminalVertex b('b');

  NonterminalVertex ba(b, a);
  NonterminalVertex bab(ba, b);
  NonterminalVertex baba(bab, a);

  FilledMatchingTable table(ba, baba);
  EXPECT_EQ(FiniteArithmeticSequence(2, 1, 1), table.matches(ba, baba));
}

TEST(FilledMatchingTable, Example3) {
  TerminalVertex a('a');
  TerminalVertex b('b');

  NonterminalVertex ba(b, a);
  NonterminalVertex ab(a, b);
  NonterminalVertex baab(ba, ab);

  FilledMatchingTable table(baab, baab);
  EXPECT_EQ(FiniteArithmeticSequence(0, 1, 1), table.matches(baab, baab));
}

TEST(LocalSearch, Trivial) {
  TerminalVertex a('a');

  FilledMatchingTable matching_table(a, a);
  PatternMatchesGenerator inspector(a, a, 0, 1, &matching_table);

  EXPECT_EQ(FiniteArithmeticSequence(0, 1, 1), inspector.next_match());
  EXPECT_EQ(FiniteArithmeticSequence(), inspector.next_match());
}

TEST(LocalSearch, LeftCut) {
  TerminalVertex a('a');
  NonterminalVertex a2(a, a);

  FilledMatchingTable matching_table(a, a2);
  PatternMatchesGenerator inspector(a, a2, 0, 1, &matching_table);

  EXPECT_EQ(FiniteArithmeticSequence(0, 1, 1), inspector.next_match());
  EXPECT_EQ(FiniteArithmeticSequence(), inspector.next_match());
}

TEST(LocalSearch, RightCut) {
  TerminalVertex a('a');
  NonterminalVertex a2(a, a);

  FilledMatchingTable matching_table(a, a2);
  PatternMatchesGenerator inspector(a, a2, 1, 1, &matching_table);

  EXPECT_EQ(FiniteArithmeticSequence(1, 1, 1), inspector.next_match());
  EXPECT_EQ(FiniteArithmeticSequence(), inspector.next_match());
}

TEST(LocalSearch, SimpleNontrivial) {
  TerminalVertex a('a');
  TerminalVertex b('b');
  NonterminalVertex ab(a, b);
  NonterminalVertex abab(ab, ab);

  FilledMatchingTable matching_table(ab, abab);
  PatternMatchesGenerator inspector(ab, abab, 0, 4, &matching_table);

  EXPECT_EQ(FiniteArithmeticSequence(0, 2, 2), inspector.next_match());
  EXPECT_EQ(FiniteArithmeticSequence(), inspector.next_match());
}

TEST(LocalSearch, SimpleNontrivialLeftCut) {
  TerminalVertex a('a');
  TerminalVertex b('b');
  NonterminalVertex ab(a, b);
  NonterminalVertex abab(ab, ab);

  FilledMatchingTable matching_table(ab, abab);
  PatternMatchesGenerator inspector(ab, abab, 1, 4, &matching_table);

  EXPECT_EQ(FiniteArithmeticSequence(2, 1, 1), inspector.next_match());
  EXPECT_EQ(FiniteArithmeticSequence(), inspector.next_match());
}

TEST(LocalSearch, SimpleNontrivialRightCut) {
  TerminalVertex a('a');
  TerminalVertex b('b');
  NonterminalVertex ab(a, b);
  NonterminalVertex abab(ab, ab);

  FilledMatchingTable matching_table(ab, abab);
  PatternMatchesGenerator inspector(ab, abab, 0, 3, &matching_table);
  EXPECT_EQ(FiniteArithmeticSequence(0, 1, 1), inspector.next_match());
  EXPECT_EQ(FiniteArithmeticSequence(), inspector.next_match());
}

TEST(LocalSearch, SimpleNontrivialSplitted) {
  TerminalVertex a('a');
  TerminalVertex b('b');
  NonterminalVertex ab(a, b);
  NonterminalVertex aba(ab, a);
  NonterminalVertex abaab(aba, ab);

  FilledMatchingTable ab_matching_table(ab, abaab);
  PatternMatchesGenerator ab_inspector(ab, abaab, 0, 5, &ab_matching_table);


  EXPECT_EQ(FiniteArithmeticSequence(0, 1, 1), ab_inspector.next_match());
  EXPECT_EQ(FiniteArithmeticSequence(3, 1, 1), ab_inspector.next_match());
  EXPECT_EQ(FiniteArithmeticSequence(), ab_inspector.next_match());

  FilledMatchingTable a_matching_table(a, abaab);
  PatternMatchesGenerator a_inspector(a, abaab, 0, 5, &a_matching_table);

  EXPECT_EQ(FiniteArithmeticSequence(0, 1, 1), a_inspector.next_match());
  EXPECT_EQ(FiniteArithmeticSequence(2, 1, 2), a_inspector.next_match());
  EXPECT_EQ(FiniteArithmeticSequence(), a_inspector.next_match());
}

Vertex get_random_slp_on_2_letters(unsigned int WORD_SIZE) {
  TerminalVertex a('a');
  TerminalVertex b('b');

  int random_word = rand() % (1 << WORD_SIZE);
  std::vector<unsigned int> random_word_split;
  for (unsigned int i = 1; i < WORD_SIZE; ++i) {
    random_word_split.push_back(i);
  }

  std::random_shuffle(random_word_split.begin(), random_word_split.end());

  std::vector<Vertex> word_presentation;
  for (unsigned int i = 0; i < WORD_SIZE; ++i) {
    word_presentation.push_back(((random_word & (1 << i)) ? b : a));
  }

  for (unsigned int split : random_word_split) {
    NonterminalVertex new_vertex(word_presentation[split - 1], word_presentation[split]);
    for (unsigned int i = split - new_vertex.left_child().length().get_ui(); i < split + new_vertex.right_child().length(); ++i) {
      word_presentation[i] = new_vertex;
    }
  }

  return word_presentation.front();
}

std::string print_tree_preorder(const Vertex& vertex) {
  std::ostringstream out;
  PreorderInspector inspector(vertex);
  while (!inspector.stopped()) {
    PrintTo(inspector.vertex(), &out);
    out << " -- ";

    ++inspector;
  }

  return out.str();
}

TEST(LocalSearch, RandomWord) {
  const unsigned int WORD_SIZE = 16;
  int REPEAT = 10000;

  srand(time(NULL));

  while (--REPEAT >= 0) {
    Vertex text = get_random_slp_on_2_letters(WORD_SIZE);

    int random_pattern_number = rand() % (2 * WORD_SIZE - 1);
    PostorderInspector pattern_getter(text);
    int i = random_pattern_number;
    while (--i >= 0) {
      ++pattern_getter;
    }

    Vertex pattern = pattern_getter.vertex();

    unsigned int left_boundary = WORD_SIZE;
    unsigned int right_boundary = 0;
    while (left_boundary >= right_boundary) {
      left_boundary = rand() % WORD_SIZE;
      right_boundary = rand() % WORD_SIZE + 1;
    }

    unsigned int text_length = right_boundary - left_boundary;

    FilledMatchingTable matching_table(pattern, text);
    PatternMatchesGenerator generator(pattern, text, left_boundary, text_length, &matching_table);

    auto pattern_tree_string = print_tree_preorder(pattern);
    auto text_tree_string = print_tree_preorder(text);

    ::std::string pattern_string(VertexWord(pattern).begin(), VertexWord(pattern).end());
    ::std::string text_string(VertexWord(text).begin(), VertexWord(text).end());

    unsigned int last_match_position = left_boundary;

    auto match = generator.next_match();
    int last_approved_match = -1;
    while (match != FiniteArithmeticSequence()) {
      unsigned int checked_matches = 0;
      unsigned int next_match_to_check = match.first().get_ui();
      while (next_match_to_check <= match.last()) {
        if (last_approved_match != next_match_to_check) {
          ASSERT_LE(next_match_to_check + pattern.length().get_ui(), right_boundary)
              << "Local search found match of " << pattern_string
              << " in " << text_string << "[" << left_boundary << ":" << right_boundary
              << "] out of boundaries. Text: " << text_tree_string
              << " pattern: " << pattern_tree_string;

          last_match_position = text_string.find(pattern_string, last_match_position);

          ASSERT_NE(std::string::npos, last_match_position)
              << "Naive algorithm can't find any more matches of " << pattern_string
              << " in " << text_string << "[" << left_boundary << ":" << right_boundary
              << "], while local search found " << next_match_to_check
              << " Text: " << text_tree_string
              << " pattern:" << pattern_tree_string;

          ASSERT_EQ(next_match_to_check, last_match_position)
              << "Naive algorithm can't find match " << next_match_to_check << " of " << pattern_string
              << " in " << text_string << "[" << left_boundary << ":" << right_boundary
              << "], found " << last_match_position << " instead. Text: " << text_tree_string
              << " pattern: " << pattern_tree_string;
          ++last_match_position;
        }
        ++checked_matches;
        last_approved_match = next_match_to_check;
        next_match_to_check += match.step().get_ui();
      }
      match = generator.next_match();
    }

    last_match_position = text_string.find(pattern_string, last_match_position);

    ASSERT_TRUE(last_match_position == std::string::npos || last_match_position + pattern.length().get_ui() > right_boundary)
      << "Naive algorithm found one more match of " << pattern_string
      << " in " << text_string << "[" << left_boundary << ":" << right_boundary
      << "] at position" << last_match_position
      << " Text: " << text_tree_string
      << " pattern: " << pattern_tree_string;
  }
}

TEST(SLPMatchingTable, Trivial) {
  TerminalVertex a('a');
  TerminalVertex a_copy('a');
  NonterminalVertex a2(a, a_copy);

  MatchingTable matching_table;

  EXPECT_EQ(FiniteArithmeticSequence(0, 1, 1), matching_table.matches(a, a));
  EXPECT_EQ(FiniteArithmeticSequence(0, 1, 1), matching_table.matches(a, a_copy));
  EXPECT_EQ(FiniteArithmeticSequence(0, 1, 2), matching_table.matches(a, a2));
}

TEST(SLPMatchingTable, Example1) {
  TerminalVertex a('a');
  TerminalVertex b('b');
  NonterminalVertex ab(a, b);

  MatchingTable matching_table;

  EXPECT_EQ(FiniteArithmeticSequence(0, 1, 1), matching_table.matches(a, ab));
  EXPECT_EQ(FiniteArithmeticSequence(1, 1, 1), matching_table.matches(b, ab));
  EXPECT_EQ(FiniteArithmeticSequence(0, 1, 1), matching_table.matches(ab, ab));
}

TEST(SLPMatchingTable, Example2) {
  TerminalVertex a('a');
  TerminalVertex b('b');
  NonterminalVertex ba(b, a);
  NonterminalVertex bb(b, b);
  NonterminalVertex bba(bb, a);

  MatchingTable matching_table;

  EXPECT_EQ(FiniteArithmeticSequence(), matching_table.matches(ba, bb));
  EXPECT_EQ(FiniteArithmeticSequence(1, 1, 1), matching_table.matches(ba, bba));
}

TEST(SLPMatchingTable, Example3) {
  TerminalVertex a('a');
  TerminalVertex b('b');

  NonterminalVertex ba(b, a);
  NonterminalVertex bab(ba, b);
  NonterminalVertex baba(bab, a);

  MatchingTable table;
  EXPECT_EQ(FiniteArithmeticSequence(2, 1, 1), table.matches(ba, baba));
}

TEST(SLPMatchingTable, Example4) {
  TerminalVertex a('a');
  TerminalVertex b('b');

  NonterminalVertex bb(b, b);
  NonterminalVertex bbb(bb, b);
  NonterminalVertex abbb(a, bbb);

  MatchingTable table;
  EXPECT_EQ(FiniteArithmeticSequence(1, 1, 1), table.matches(bb, abbb));
}

TEST(SLPMatchingTable, Example5) {
  TerminalVertex a('a');
  TerminalVertex b('b');

  NonterminalVertex aa(a, a);
  NonterminalVertex aab(aa, b);
  NonterminalVertex baab(b, aab);

  MatchingTable table;
  EXPECT_EQ(FiniteArithmeticSequence(0, 1, 1), table.matches(aa, aa));
  EXPECT_EQ(FiniteArithmeticSequence(0, 1, 1), table.matches(aa, aab));
  EXPECT_EQ(FiniteArithmeticSequence(1, 1, 1), table.matches(aa, baab));
  EXPECT_EQ(FiniteArithmeticSequence(1, 1, 1), table.matches(aab, baab));
}

TEST(SLPMatchingTable, Example6) {
  TerminalVertex a('a');

  NonterminalVertex aa(a, a);
  NonterminalVertex aaa(aa, a);
  NonterminalVertex a5(aa, aaa);
  NonterminalVertex a6(a5, a);

  MatchingTable table;
  EXPECT_EQ(FiniteArithmeticSequence(0, 1, 1), table.matches(aa, aa));
  EXPECT_EQ(FiniteArithmeticSequence(0, 1, 2), table.matches(aa, aaa));
  EXPECT_EQ(FiniteArithmeticSequence(0, 1, 3), table.matches(aa, a5));
  EXPECT_EQ(FiniteArithmeticSequence(3, 1, 2), table.matches(aa, a6));
  EXPECT_EQ(FiniteArithmeticSequence(0, 1, 1), table.matches(aaa, aaa));
  EXPECT_EQ(FiniteArithmeticSequence(0, 1, 3), table.matches(aaa, a5));
  EXPECT_EQ(FiniteArithmeticSequence(2, 1, 2), table.matches(aaa, a6));
  EXPECT_EQ(FiniteArithmeticSequence(0, 1, 1), table.matches(a5, a5));
  EXPECT_EQ(FiniteArithmeticSequence(0, 1, 2), table.matches(a5, a6));
}

TEST(SLPMatchingTable, Example7) {
  TerminalVertex a('a');
  TerminalVertex b('b');

  NonterminalVertex v1(a, b);
  NonterminalVertex v2(b, b);
  NonterminalVertex v3(b, b);
  NonterminalVertex v4(b, v2);
  NonterminalVertex v5(v1, b);
  NonterminalVertex v6(v3, v4);
  NonterminalVertex v7(v5, v6);

  MatchingTable table;
  EXPECT_EQ(FiniteArithmeticSequence(), table.matches(v2, v1));
  EXPECT_EQ(FiniteArithmeticSequence(1, 1, 1), table.matches(v2, v5));
  EXPECT_EQ(FiniteArithmeticSequence(), table.matches(v4, v5));
  EXPECT_EQ(FiniteArithmeticSequence(1, 1, 3), table.matches(v2, v7));
  EXPECT_EQ(FiniteArithmeticSequence(0, 1, 1), table.matches(v2, v3));
  EXPECT_EQ(FiniteArithmeticSequence(0, 1, 3), table.matches(v2, v6));
  EXPECT_EQ(FiniteArithmeticSequence(1, 1, 3), table.matches(v4, v7));
  EXPECT_EQ(FiniteArithmeticSequence(0, 1, 2), table.matches(v2, v4));
  EXPECT_EQ(FiniteArithmeticSequence(0, 1, 3), table.matches(v4, v6));
  EXPECT_EQ(FiniteArithmeticSequence(), table.matches(v3, v1));
  EXPECT_EQ(FiniteArithmeticSequence(1, 1, 1), table.matches(v3, v5));
  EXPECT_EQ(FiniteArithmeticSequence(1, 1, 3), table.matches(v3, v7));
  EXPECT_EQ(FiniteArithmeticSequence(0, 1, 3), table.matches(v3, v6));
  EXPECT_EQ(FiniteArithmeticSequence(1, 1, 3), table.matches(v6, v7));
}

TEST(SLPMatchingTable, InversedTest) {
  TerminalVertex a('a');
  TerminalVertex b('b');

  NonterminalVertex ab(a, b);
  NonterminalVertex a_1b(a.negate(), b);
  NonterminalVertex b_1a_1b(b.negate(), a_1b);

  MatchingTable matching_table;

  EXPECT_EQ(FiniteArithmeticSequence(1, 1, 1), matching_table.matches(a, b_1a_1b.negate()));
  EXPECT_EQ(FiniteArithmeticSequence(2, 1, 1), matching_table.matches(b, b_1a_1b.negate()));
  EXPECT_EQ(FiniteArithmeticSequence(0, 1, 1), matching_table.matches(b.negate(), a_1b.negate()));
  EXPECT_EQ(FiniteArithmeticSequence(1, 1, 1), matching_table.matches(ab, b_1a_1b.negate()));
  EXPECT_EQ(FiniteArithmeticSequence(0, 1, 1), matching_table.matches(ab.negate(), b_1a_1b));
}

TEST(SLPMatchingTable, StressTest) {
  const unsigned int WORD_SIZE = 10;
  int repeat = 5000;

  srand(time(NULL));

  while (--repeat >= 0) {
    Vertex slp = get_random_slp_on_2_letters(WORD_SIZE);

    MatchingTable real_matching_table;
    FilledMatchingTable computed_matching_table(slp, slp);

    PostorderInspector text(slp);
    unsigned int text_vertex_number = 0;

    while (!text.stopped()) {
      ++text_vertex_number;

      PostorderInspector pattern(slp);
      unsigned int pattern_vertex_number = 0;
      while (!pattern.stopped()) {
        ++pattern_vertex_number;
        ASSERT_EQ(computed_matching_table.matches(pattern.vertex(), text.vertex()), real_matching_table.matches(pattern.vertex(), text.vertex()))
          << "pattern " << VertexWord(pattern.vertex()) << ::std::endl << print_tree_preorder(pattern.vertex()) << ::std::endl << ::std::endl
          << "text " << VertexWord(text.vertex()) << ::std::endl << print_tree_preorder(text.vertex()) << ::std::endl << ::std::endl;
        pattern.next();
      }
      text.next();
    }
  }
}

}
}
}


