/**
 * \file slp_common_prefix.cpp
 * \brief Tests for slp_common_prefix.h
 *
 *  Created on: Jan 21, 2013
 *      Author: dpantele
 */

#include <memory>
#include <functional>
#include <utility>
#include <string>
#include <algorithm>

#include "gtest/gtest.h"
#include "slp.h"

namespace crag {
namespace slp {
namespace {

typedef TerminalVertexTemplate<char> TerminalVertex;

TEST(CommonPrefix, Trivial) {
  TerminalVertex a('a');
  NonterminalVertex aa(a, a);

  EXPECT_EQ(1, longest_common_prefix(a, a));
  EXPECT_EQ(2, longest_common_prefix(aa, aa));
  EXPECT_EQ(1, longest_common_prefix(a, aa));
}

TEST(CommonPrefix, Example1) {
  TerminalVertex a('a');
  NonterminalVertex aa(a, a);
  NonterminalVertex a3left(aa, a);
  NonterminalVertex a3right(a, aa);

  EXPECT_EQ(1, longest_common_prefix(a, a3left));
  EXPECT_EQ(1, longest_common_prefix(a, a3right));
  EXPECT_EQ(2, longest_common_prefix(aa, a3left));
  EXPECT_EQ(2, longest_common_prefix(aa, a3right));
  EXPECT_EQ(3, longest_common_prefix(a3left, a3right));
}

TEST(CommonPrefix, Example2) {
  TerminalVertex a('a');
  TerminalVertex b('b');
  NonterminalVertex aa(a, a);
  NonterminalVertex ab(a, b);
  MatchingTable matching_table;
  matching_table.matches(aa, ab);

  EXPECT_EQ(1, longest_common_prefix(aa, ab));
  EXPECT_EQ(1, longest_common_prefix(aa, ab, &matching_table));
}

TEST(CommonPrefix, Example3) {
  TerminalVertex a('a');
  TerminalVertex b('b');
  NonterminalVertex ab(a, b);
  NonterminalVertex b1ab(ab.negate(), b);
  NonterminalVertex bb1ab(b, b1ab);

  MatchingTable matching_table;

  EXPECT_EQ(0, get_cancellation_length(a, &matching_table));
  EXPECT_EQ(0, get_cancellation_length(b, &matching_table));
  EXPECT_EQ(0, get_cancellation_length(ab, &matching_table));
  EXPECT_EQ(0, get_cancellation_length(b1ab, &matching_table));
//  std::cout << matching_table.matches(b, ab) << std::endl;
//  std::cout << matching_table.matches(b.negate(), ab.negate()) << std::endl;
  EXPECT_EQ(1, get_cancellation_length(bb1ab, &matching_table));
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

Vertex get_prefix(const Vertex& slp, LongInteger prefix_length) {
  assert(slp.height() > 0);
  if (slp.height() == 1) {
    assert(prefix_length == 1);
    return slp;
  }
  if (prefix_length <= slp.split_point()) {
    return get_prefix(slp.left_child(), prefix_length);
  }
  return NonterminalVertex(slp.left_child(), get_prefix(slp.right_child(), prefix_length - slp.split_point()));
}

TEST(GetPrefix, StressTest) {
  const size_t size = 10;
  int repeat = 10000;
  while (--repeat >= 0) {
    Vertex slp = get_random_slp_on_2_letters(size);
    std::string slp_string(VertexWord<char>(slp).begin(), VertexWord<char>(slp).end());
    for (size_t prefix_length = 1; prefix_length < size; ++prefix_length) {
      Vertex prefix = get_prefix(slp, prefix_length);
      ASSERT_EQ(prefix_length, prefix.length());
      std::string prefix_string(VertexWord<char>(prefix).begin(), VertexWord<char>(prefix).end());
      ASSERT_EQ(slp_string.substr(0, prefix_length), prefix_string);
    }
  }
}

TEST(CommonPrefix, StressTest) {
  const size_t size = 10;
  int repeat = 1000;
  srand(time(0));
  while (--repeat >= 0) {
    Vertex slp = get_random_slp_on_2_letters(size);
    for (size_t prefix_length = 1; prefix_length < size; ++prefix_length) {
      Vertex prefix = get_prefix(slp, prefix_length);
      ASSERT_EQ(prefix_length, longest_common_prefix(prefix, slp));
    }
  }
}

TEST(CommonPrefix, StressTestWithCommonTable) {
  const size_t size = 10;
  int repeat = 1000;
  srand(time(0));
  MatchingTable common_table;
  while (--repeat >= 0) {
    Vertex slp = get_random_slp_on_2_letters(size);
    for (size_t prefix_length = 1; prefix_length < size; ++prefix_length) {
      Vertex prefix = get_prefix(slp, prefix_length);
      ASSERT_EQ(prefix_length, longest_common_prefix(prefix, slp, &common_table));
    }
  }
}

} //namespace
} //slp
} //crag


