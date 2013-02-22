/*
 * SLP_vertex_word.cpp
 *
 *  Created on: Feb 21, 2013
 *      Author: dpantele
 */

#include "gtest/gtest.h"
#include "slp.h"

#include <memory>
#include <functional>
#include <utility>
#include <algorithm>
#include <string>

namespace crag {
namespace slp {
namespace {
  typedef VertexWord<char> VertexCharWord;
  typedef TerminalVertexTemplate<char> TerminalVertex;

  TEST(VertexWord, size) {
    TerminalVertex a('a');
    TerminalVertex b('a');
    NonterminalVertex ab(a, b);

    EXPECT_EQ(0, VertexCharWord().size());
    EXPECT_EQ(1, VertexCharWord(a).size());
    EXPECT_EQ(2, VertexCharWord(ab).size());
  }

  TEST(VertexWord, IndexLeftTree) {
    TerminalVertex a('a');
    TerminalVertex b('b');
    TerminalVertex c('c');

    NonterminalVertex ab(a, b);
    NonterminalVertex abc(ab, c);

    VertexCharWord w(abc);

    ASSERT_EQ(3, w.size());
    EXPECT_EQ('a', w[0]);
    EXPECT_EQ('b', w[1]);
    EXPECT_EQ('c', w[2]);
  }

  TEST(VertexWord, IndexRightTree) {
    TerminalVertex a('a');
    TerminalVertex b('b');
    TerminalVertex c('c');

    NonterminalVertex bc(b, c);
    NonterminalVertex abc(a, bc);

    VertexCharWord w(abc);

    ASSERT_EQ(3, w.size());
    EXPECT_EQ('a', w[0]);
    EXPECT_EQ('b', w[1]);
    EXPECT_EQ('c', w[2]);
  }

  TEST(VertexWord, IndexCrossedTree) {
    TerminalVertex a('a');
    TerminalVertex b('b');

    NonterminalVertex ab(a, b);
    NonterminalVertex ba_1(b, a.negate());

    NonterminalVertex abab_1(ab, ba_1.negate());

    VertexCharWord w(abab_1);

    ASSERT_EQ(w.size(), 4);
    EXPECT_EQ('a', w[0]);
    EXPECT_EQ('b', w[1]);
    EXPECT_EQ('a', w[2]);
    EXPECT_EQ(-'b', w[3]);
  }

  TEST(VertexWord, Iterator) {
    TerminalVertex a('a');
    TerminalVertex b('b');
    TerminalVertex c('c');

    NonterminalVertex ab(a, b);
    NonterminalVertex abc(ab, c);

    VertexCharWord w(abc);

    EXPECT_EQ("abc", ::std::string(w.begin(), w.end()));
  }


  TEST(VertexWord, StressTest) {
    const unsigned int word_size = 6;

    unsigned int current_word = 0;
    while (current_word < (1 << word_size)) {
      std::vector<unsigned int> current_word_split;
      std::string current_word_string;

      for (size_t i = 0; i < word_size; ++i) {
        current_word_string.push_back((current_word & (1 << i)) ? 'b' : 'a');
      }

      for (unsigned int i = 1; i < word_size; ++i) {
        current_word_split.push_back(i);
      }

      do {
        std::vector<Vertex> word_presentation;
        for (char letter : current_word_string) {
          word_presentation.push_back(TerminalVertex(letter));
        }

        for (unsigned int split : current_word_split) {
          NonterminalVertex new_vertex(word_presentation[split - 1], word_presentation[split]);
          for (unsigned int i = split - new_vertex.left_child().length().get_ui(); i < split + new_vertex.right_child().length(); ++i) {
            word_presentation[i] = new_vertex;
          }
        }

        unsigned int position = 0;
        const Vertex& current_slp = word_presentation.front();
        ASSERT_EQ(current_word_split.back(), current_slp.left_child().length().get_ui());

        VertexCharWord word(current_slp);
        ASSERT_EQ(current_word_string, std::string(word.begin(), word.end()));

      } while(std::next_permutation(current_word_split.begin(), current_word_split.end()));
      ++current_word;
    }
  }
}
}
}


