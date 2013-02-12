/*
 * SLPSet_produced_word.cpp
 *
 *  Created on: Jan 23, 2013
 *      Author: dpantele
 */

/*
 * SLPSet_common.cpp
 *
 *  Created on: Jan 21, 2013
 *      Author: dpantele
 */

#include "gtest/gtest.h"
#include "SLPSet.h"

#include <memory>
#include <functional>
#include <utility>
#include <algorithm>

namespace crag {
namespace {
  TEST(SLPProducedWord, size) {
    EXPECT_EQ(SLPProducedWord().size(), 0);
    SLPVertex a = SLPVertex::terminal_vertex(1);
    SLPVertex b = SLPVertex::terminal_vertex(2);
    SLPVertex ab = SLPVertex::concatenate(a, b);

    EXPECT_EQ(SLPProducedWord(a).size(), 1);
    EXPECT_EQ(SLPProducedWord(ab).size(), 2);
  }

  TEST(SLPProducedWord, IndexLeftTree) {
    SLPVertex a = SLPVertex::terminal_vertex(1);
    SLPVertex b = SLPVertex::terminal_vertex(2);
    SLPVertex c = SLPVertex::terminal_vertex(3);

    SLPVertex ab = SLPVertex::concatenate(a, b);
    SLPVertex abc = SLPVertex::concatenate(ab, c);

    SLPProducedWord w(abc);

    ASSERT_EQ(w.size(), 3);
    EXPECT_EQ(w[0], a);
    EXPECT_EQ(w[1], b);
    EXPECT_EQ(w[2], c);
  }

  TEST(SLPProducedWord, IndexRightTree) {
    SLPVertex a = SLPVertex::terminal_vertex(1);
    SLPVertex b = SLPVertex::terminal_vertex(2);
    SLPVertex c = SLPVertex::terminal_vertex(3);

    SLPVertex bc = SLPVertex::concatenate(b, c);
    SLPVertex abc = SLPVertex::concatenate(a, bc);

    SLPProducedWord w(abc);

    ASSERT_EQ(w.size(), 3);
    EXPECT_EQ(w[0], a);
    EXPECT_EQ(w[1], b);
    EXPECT_EQ(w[2], c);
  }

  TEST(SLPProducedWord, IndexCrossedTree) {
    SLPVertex a = SLPVertex::terminal_vertex(1);
    SLPVertex b = SLPVertex::terminal_vertex(2);

    SLPVertex ab = SLPVertex::concatenate(a, b);
    SLPVertex ba_1 = SLPVertex::concatenate(b, a.negate());//ba^-1

    SLPVertex abab_1 = SLPVertex::concatenate(ab, ba_1.negate());

    SLPProducedWord w(abab_1);

    ASSERT_EQ(w.size(), 4);
    EXPECT_EQ(w[0], a);
    EXPECT_EQ(w[1], b);
    EXPECT_EQ(w[2], a);
    EXPECT_EQ(w[3], b.negate());
  }

  TEST(SLPProducedWord, Iterator) {
    SLPVertex a = SLPVertex::terminal_vertex(1);
    SLPVertex b = SLPVertex::terminal_vertex(2);
    SLPVertex c = SLPVertex::terminal_vertex(3);

    SLPVertex ab = SLPVertex::concatenate(a, b);
    SLPVertex abc = SLPVertex::concatenate(ab, c);

    SLPProducedWord w(abc);

    auto iterator = w.begin();
    EXPECT_NE(iterator, w.end());
    EXPECT_EQ(*iterator, a);
    ++iterator;
    EXPECT_NE(iterator, w.end());
    EXPECT_EQ(iterator->terminal_symbol(), b.terminal_symbol());
    ++iterator;
    EXPECT_NE(iterator, w.end());
    EXPECT_EQ(*iterator, c);
    ++iterator;

    EXPECT_EQ(iterator, w.end());
  }


  TEST(SLPProducedWord, StressTest) {
    const unsigned int word_size = 6;

    SLPVertex a = SLPVertex::terminal_vertex(1);
    SLPVertex b = SLPVertex::terminal_vertex(2);

    unsigned int current_word = 0;
    while (current_word < (1 << word_size)) {
      std::vector<unsigned int> current_word_split;
      for (unsigned int i = 1; i < word_size; ++i) {
        current_word_split.push_back(i);
      }

      do {
        std::vector<SLPVertex> word_presentation;
        for (unsigned int i = 0; i < word_size; ++i) {
          word_presentation.push_back(((current_word & (1 << i)) ? b : a));
        }

        for (unsigned int split : current_word_split) {
          SLPVertex new_vertex = SLPVertex::concatenate(word_presentation[split - 1], word_presentation[split]);
          for (unsigned int i = split - new_vertex.left_child().length().get_ui(); i < split + new_vertex.right_child().length(); ++i) {
            word_presentation[i] = new_vertex;
          }
        }

        unsigned int position = 0;
        EXPECT_EQ(current_word_split.back(), word_presentation.front().left_child().length());
        for (auto symbol = SLPProducedWordIterator(word_presentation.front()); symbol != SLPProducedWordIterator(); ++position, ++symbol) {
          EXPECT_EQ(((current_word & (1 << position)) ? b : a), *symbol);
        }

      } while(std::next_permutation(current_word_split.begin(), current_word_split.end()));
      ++current_word;
    }
  }
}
}

