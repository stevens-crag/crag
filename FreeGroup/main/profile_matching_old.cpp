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

#include "SLPSet.h"

using crag::SLPVertex;
using crag::SLPMatchingTable;
using crag::SLPPostorderInspector;

SLPVertex get_random_slp_on_2_letters(unsigned int WORD_SIZE) {
  SLPVertex a = SLPVertex::terminal_vertex(1);
  SLPVertex b = SLPVertex::terminal_vertex(2);

  int random_word = rand() % (1 << WORD_SIZE);
  std::vector<unsigned int> random_word_split;
  for (unsigned int i = 1; i < WORD_SIZE; ++i) {
    random_word_split.push_back(i);
  }

  std::random_shuffle(random_word_split.begin(), random_word_split.end());

  std::vector<SLPVertex> word_presentation;
  for (unsigned int i = 0; i < WORD_SIZE; ++i) {
    word_presentation.push_back(((random_word & (1 << i)) ? b : a));
  }

  for (unsigned int split : random_word_split) {
    SLPVertex new_vertex = SLPVertex::concatenate(word_presentation[split - 1], word_presentation[split]);
    for (unsigned int i = split - new_vertex.left_child().length().get_ui(); i < split + new_vertex.right_child().length().get_ui(); ++i) {
      word_presentation[i] = new_vertex;
    }
  }

  return word_presentation.front();
}

int main() {
  const unsigned int WORD_SIZE = 16;
  int REPEAT = 10000;

  srand(10010);

  auto begin = std::chrono::high_resolution_clock::now();

  while (--REPEAT >= 0) {
    SLPVertex text = get_random_slp_on_2_letters(WORD_SIZE);

    int random_pattern_number = rand() % (2 * WORD_SIZE - 1);
    SLPPostorderInspector pattern_getter(text);
    int i = random_pattern_number;
    while (--i >= 0) {
      pattern_getter.go_to_next_vertex();
    }

    SLPVertex pattern = pattern_getter.current_vertex();
    auto result = SLPMatchingTable().matches(pattern, text);
  }

  auto end = std::chrono::system_clock::now();
  std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << std::endl;
}


