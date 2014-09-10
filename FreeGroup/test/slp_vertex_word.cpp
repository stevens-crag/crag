/*
 * SLP_vertex_word.cpp
 *
 *  Created on: Feb 21, 2013
 *      Author: dpantele
 */

#include <memory>
#include <functional>
#include <utility>
#include <algorithm>
#include <string>

#include "gtest/gtest.h"
#include "slp_vertex_word.h"

namespace crag {
namespace slp {
namespace {

CONSTEXPR_OR_CONST TerminalSymbol terminal_a = TerminalSymbol{} + 1;
CONSTEXPR_OR_CONST TerminalSymbol terminal_b = TerminalSymbol{} + 2;
CONSTEXPR_OR_CONST TerminalSymbol terminal_c = TerminalSymbol{} + 3;

  TEST(VertexWord, size) {
    TerminalVertex a(terminal_a);
    TerminalVertex b(terminal_a);
    NonterminalVertex ab(a, b);

    EXPECT_EQ(0, VertexWord().size());
    EXPECT_EQ(1, VertexWord(a).size());
    EXPECT_EQ(2, VertexWord(ab).size());
  }

  TEST(VertexWord, IndexLeftTree) {
    TerminalVertex a(terminal_a);
    TerminalVertex b(terminal_b);
    TerminalVertex c(terminal_c);

    NonterminalVertex ab(a, b);
    NonterminalVertex abc(ab, c);

    VertexWord w(abc);

    ASSERT_EQ(3, w.size());
    EXPECT_EQ(terminal_a, w[0]);
    EXPECT_EQ(terminal_b, w[1]);
    EXPECT_EQ(terminal_c, w[2]);
  }

  TEST(VertexWord, IndexRightTree) {
    TerminalVertex a(terminal_a);
    TerminalVertex b(terminal_b);
    TerminalVertex c(terminal_c);

    NonterminalVertex bc(b, c);
    NonterminalVertex abc(a, bc);

    VertexWord w(abc);

    ASSERT_EQ(3, w.size());
    EXPECT_EQ(terminal_a, w[0]);
    EXPECT_EQ(terminal_b, w[1]);
    EXPECT_EQ(terminal_c, w[2]);
  }

  TEST(VertexWord, IndexCrossedTree) {
    TerminalVertex a(terminal_a);
    TerminalVertex b(terminal_b);

    NonterminalVertex ab(a, b);
    NonterminalVertex ba_1(b, a.negate());

    NonterminalVertex abab_1(ab, ba_1.negate());

    VertexWord w(abab_1);

    ASSERT_EQ(w.size(), 4);
    EXPECT_EQ(terminal_a, w[0]);
    EXPECT_EQ(terminal_b, w[1]);
    EXPECT_EQ(terminal_a, w[2]);
    EXPECT_EQ(-terminal_b, w[3]);
  }

  TEST(VertexWord, Iterator) {
    TerminalVertex a(terminal_a);
    TerminalVertex b(terminal_b);
    TerminalVertex c(terminal_c);

    NonterminalVertex ab(a, b);
    NonterminalVertex abc(ab, c);

    VertexWord w(abc);

    EXPECT_EQ(::std::string({static_cast<char>(terminal_a), static_cast<char>(terminal_b), static_cast<char>(terminal_c)}), ::std::string(w.begin(), w.end()));
    EXPECT_EQ(terminal_a, *(w.begin() + 0));
    EXPECT_EQ(terminal_b, *(w.begin() + 1));
    EXPECT_EQ(terminal_c, *(w.begin() + 2));
  }


  TEST(VertexWord, StressTest) {
    const unsigned int word_size = 6;
    TerminalVertex a(terminal_a);
    TerminalVertex b(terminal_b);

    unsigned int current_word = 0;
    while (current_word < (1 << word_size)) {
      std::vector<unsigned int> current_word_split;
      std::string current_word_string;

      for (size_t i = 0; i < word_size; ++i) {
        current_word_string.push_back((current_word & (1 << i)) ? terminal_b : terminal_a);
      }

      for (unsigned int i = 1; i < word_size; ++i) {
        current_word_split.push_back(i);
      }

      do {
        std::vector<Vertex> word_presentation;
        for (char letter : current_word_string) {
          word_presentation.push_back(letter == terminal_a ? a : b);
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

        VertexWord word(current_slp);
        ASSERT_EQ(current_word_string, std::string(word.begin(), word.end()));

        for (size_t i = 0; i < current_word_string.size(); ++i) {
          ASSERT_EQ(current_word_string[i], *(word.begin() + i)) <<
              "word.begin() + " << i << " does not refer to " << current_word_string[i];
        }

      } while(std::next_permutation(current_word_split.begin(), current_word_split.end()));
      ++current_word;
    }
  }
}
}
}


