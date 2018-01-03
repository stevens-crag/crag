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
#include <vector>
#include <iterator>
#include <sstream>

#include "gtest/gtest.h"
#include "slp_reduce.h"
#include "EndomorphismSLP.h"
#include <iostream>

namespace crag {
namespace slp {
namespace {

CONSTEXPR_OR_CONST TerminalSymbol terminal_a = TerminalSymbol{} + 1;
CONSTEXPR_OR_CONST TerminalSymbol terminal_b = TerminalSymbol{} + 2;

Vertex get_random_slp_on_2_letters(unsigned int WORD_SIZE) {
  TerminalVertex a(terminal_a);
  TerminalVertex b(terminal_b);

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
//    std::string slp_string(VertexWord(slp).begin(), VertexWord(slp).end());
//    for (size_t start = 1; start < size; ++start) {
//      for (size_t end = start + 1; end <= size; ++end) {
//        Vertex sub_slp = get_sub_slp(slp, start, end);
//        std::string subslp_string(VertexWord(sub_slp).begin(), VertexWord(sub_slp).end());
//        ASSERT_EQ(slp_string.substr(start, end - start), subslp_string) << slp_string << '[' << start << ':' <<end << ']';
//      }
//    }
//  }
//}

TEST(SubSLP, StressTest) {
  const size_t size = 10;
  int repeat = 1000;
  srand(time(0));
  while (--repeat >= 0) {
    Vertex slp = get_random_slp_on_2_letters(size);
    std::string slp_string(VertexWord(slp).begin(), VertexWord(slp).end());
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
    std::string subslp_string(VertexWord(sub_slp).begin(), VertexWord(sub_slp).end());
    ASSERT_EQ(slp_string.substr(start, end - start), subslp_string) << slp_string << '[' << start << ':' <<end << ']';

    PreorderInspector inspector(sub_slp);
    while (!inspector.stopped()) {
      EXPECT_TRUE(inspector.vertex().height() == 1 || (inspector.vertex().left_child() && inspector.vertex().right_child()))
          << "Vertex " << ::testing::PrintToString(inspector.vertex()) << " has one empty child";
      ++inspector;
    }
  }
}

TEST(Reduce, Simple1) {
  TerminalVertex a(terminal_a);
  TerminalVertex a_(-terminal_a);
  NonterminalVertex mNull(a, a_);

  EXPECT_EQ(Vertex(), reduce(mNull));
}

TEST(Reduce, Simple2) {
  TerminalVertex a(terminal_a);
  TerminalVertex b(terminal_b);
  NonterminalVertex ab(a, b);
  NonterminalVertex b_1a_1(b.negate(), a.negate());
  NonterminalVertex mNull(ab, b_1a_1);

  EXPECT_EQ(Vertex(), reduce(mNull));
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

std::string print_tree_preorder_single(const Vertex& vertex) {
  std::ostringstream out;
  std::unordered_map<slp::Vertex, bool> mapping;
  auto acceptor = [&] (const inspector::InspectorTask& task) {
    return mapping.find(task.vertex) == mapping.end();//do not accept if vertex is mapped already
  };

  Inspector<inspector::Preorder, decltype(acceptor)> inspector(vertex, acceptor);
  while (!inspector.stopped()) {
    PrintTo(inspector.vertex(), &out);
    out << " << ";
    if (mapping.count(inspector.vertex().left_child()) || mapping.count(inspector.vertex().left_child().negate())) {
      out << "(p) ";
    }
    PrintTo(inspector.vertex().left_child(), &out);
    out << " >> ";
    if (mapping.count(inspector.vertex().right_child()) || mapping.count(inspector.vertex().right_child().negate())) {
      out << "(p) ";
    }
    PrintTo(inspector.vertex().right_child(), &out);
    out << std::endl;
    mapping.insert(std::make_pair(inspector.vertex(), true));

    ++inspector;
  }

  return out.str();
}


TEST(Reduce, Example1) {
  TerminalVertex v1(1);
  TerminalVertex v2(2);
  NonterminalVertex v9(v1, v2);
  Vertex v12 = NonterminalVertex(v2, v9).negate();
  NonterminalVertex v13(v2, v12);
  auto reduced = reduce(v13);
  int previous = 0;
  for (auto symbol : VertexWord(reduced)) {
    EXPECT_NE(-previous, symbol) << print_tree_preorder(reduced);
    previous = symbol;
  }
}
//TerminalVertex t1(1);
//TerminalVertex t2(2);
//TerminalVertex t3(3);
//NonterminalVertex v1(t1, t2);
//NonterminalVertex v2(v1, t3);
//NonterminalVertex v3(t1, v2);
//NonterminalVertex v4(v3, v2);
//NonterminalVertex v5(v4, t3);
//NonterminalVertex v6(v5, v4);
//NonterminalVertex v7(v6.negate(), v5);
//NonterminalVertex v8(v2, t3);
//NonterminalVertex v9(v8, t3);
//NonterminalVertex v10(v7, v9);
//
//std::vector<int> correct = {-3, -2, -1, -3, -2, -1, 2, 3, 3, 3};

TEST(Reduce, Example2) {
  TerminalVertex t1(1);
  TerminalVertex t2(2);
  TerminalVertex t3(3);
  NonterminalVertex v1(t1, t2);
  NonterminalVertex v2(v1, t3);
  NonterminalVertex v3(t1, v2);
  NonterminalVertex v4(v3, v2);
  NonterminalVertex v5(v4, t3);
  NonterminalVertex v6(v5, v4);
  NonterminalVertex v7(v6.negate(), v5);
  NonterminalVertex v8(v7, v2);

  std::vector<int> correct = {-3, -2, -1, -3, -2, -1, 2, 3};

  auto reduced = reduce(v8);
  ASSERT_EQ(correct.size(), reduced.length());
  auto correct_symbol = correct.begin();
  for (auto symbol : VertexWord(reduced)) {
    ASSERT_EQ(*correct_symbol, symbol);
    ++correct_symbol;
  }
}


TEST(Reduce, StressTest) {
  const size_t REPEAT = 10000;
  const int RANK = 3;
  const size_t ENDOMORPHISMS_NUMBER = 20;
  size_t seed = 0;
  while (++seed <= REPEAT) {
    UniformAutomorphismSLPGenerator<> generator(RANK, seed);
    auto image = EndomorphismSLP::composition(ENDOMORPHISMS_NUMBER, generator).image(1);

    Vertex reduced = reduce(image);

    std::vector<int> reduced_image;
    for (auto symbol : VertexWord(reduced)) {
      if (!reduced_image.empty() && symbol == -reduced_image.back()) {
        reduced_image.pop_back();
      } else {
        reduced_image.push_back(symbol);
      }
    }
    std::ostringstream reduced_image_string;
    std::copy(reduced_image.begin(), reduced_image.end(), std::ostream_iterator<int>(reduced_image_string, ""));
    auto correct_symbol = reduced_image.begin();
    for (auto symbol : VertexWord(reduced)) {
      ASSERT_EQ(*correct_symbol, symbol) << seed << std::endl
          << print_tree_preorder_single(image) << std::endl
          << print_tree_preorder(reduced) << std::endl
          << VertexWord(image) << std::endl
          << VertexWord(reduced) << std::endl
          << reduced_image_string.str() << std::endl;
      ++correct_symbol;
    }
  }
}

//TEST(Reduce, PerformanceTest) {
//  int REPEAT = 10;
//  const size_t RANK = 3;
//  const size_t ENDOMORPHISMS_NUMBER = 100;
//  size_t seed = 112233;///time(0);//
//  UniformAutomorphismSLPGenerator<> generator(RANK, seed);
//  while (--REPEAT >= 0) {
//    auto image = EndomorphismSLP::composition(ENDOMORPHISMS_NUMBER, generator).slp(1);
//
//    Vertex reduced = reduce(image);
//    std::cout << image.length() << std::endl;
//    std::cout << reduced.length() << std::endl;
//  }
//}


} //namespace
} //slp
} //crag


