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

//TEST(SubSLP, StressTestNonRecursive) {
//  const size_t size = 2;
//  int repeat = 1;
//  while (--repeat >= 0) {
//    Vertex slp = get_random_slp_on_2_letters(size);
//    std::string slp_string(VertexWord<char>(slp).begin(), VertexWord<char>(slp).end());
//    for (size_t start = 1; start < size; ++start) {
//      for (size_t end = start + 1; end <= size; ++end) {
//        Vertex sub_slp = get_sub_slp(slp, start, end);
//        std::string subslp_string(VertexWord<char>(sub_slp).begin(), VertexWord<char>(sub_slp).end());
//        ASSERT_EQ(slp_string.substr(start, end - start), subslp_string) << slp_string << '[' << start << ':' <<end << ']';
//      }
//    }
//  }
//}

TEST(SubSLP, StressTest) {
  const size_t size = 500;
  int repeat = 1000;
  srand(time(0));
  while (--repeat >= 0) {
    Vertex slp = get_random_slp_on_2_letters(size);
    std::string slp_string(VertexWord<char>(slp).begin(), VertexWord<char>(slp).end());
    size_t start = 0;
    size_t end = 0;
    while (start == end) {
      start = rand() % size;
      end = rand() % size + 1;
      if (start > end) {
        std::swap(start, end);
      }
    }
    Vertex sub_slp = get_sub_slp(slp, start, end);
    std::string subslp_string(VertexWord<char>(sub_slp).begin(), VertexWord<char>(sub_slp).end());
    ASSERT_EQ(slp_string.substr(start, end - start), subslp_string) << slp_string << '[' << start << ':' <<end << ']';

    PreorderInspector inspector(sub_slp);
    while (!inspector.stopped()) {
      EXPECT_TRUE(inspector.vertex().height() == 1 || (inspector.vertex().left_child() && inspector.vertex().right_child()))
          << "Vertex " << ::testing::PrintToString(inspector.vertex()) << " has one empty child";
      ++inspector;
    }
  }
}


} //namespace
} //slp
} //crag


