/*
 * profile_matching_new.cpp
 *
 *  Created on: Apr 2, 2013
 *      Author: dpantele
 */

#include <memory>
#include <functional>
#include <utility>
#include <vector>
#include <tuple>
#include <algorithm>
#include <chrono>

#include "slp_inspector.h"
#include "slp_pattern_matching.h"
#include "slp_vertex_word.h"
#include "gmp_boost_pool_allocator.h"

using crag::slp::Vertex;
using crag::slp::PreorderInspector;
typedef crag::slp::TerminalVertexTemplate<char> TerminalVertex;
using crag::slp::NonterminalVertex;
using crag::slp::MatchingTable;

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
    for (unsigned int i = split - new_vertex.left_child().length().get_ui(); i < split + new_vertex.right_child().length().get_ui(); ++i) {
      word_presentation[i] = new_vertex;
    }
  }

  return word_presentation.front();
}

int main() {
  gmp_pool_setup();
  const unsigned int WORD_SIZE = 16;
  int REPEAT = 10000;

  srand(10010);

  auto begin = std::chrono::high_resolution_clock::now();

  while (--REPEAT >= 0) {
    Vertex text = get_random_slp_on_2_letters(WORD_SIZE);

    int random_pattern_number = rand() % (2 * WORD_SIZE - 1);
    PreorderInspector pattern_getter(text);
    int i = random_pattern_number;
    while (--i >= 0) {
      ++pattern_getter;
    }

    Vertex pattern = pattern_getter.vertex();
    auto result = MatchingTable().matches(pattern, text);
  }

  auto end = std::chrono::system_clock::now();
  std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << std::endl;
}


