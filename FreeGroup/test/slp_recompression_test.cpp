/**
 * \file slp_recompression.cpp
 * \brief Tests for slp_recompression.h
 */

#include <iostream>

#include "gtest/gtest.h"
#include "slp_recompression.h"
#include "slp_vertex_word.h"
#include "slp_vertex_hash.h"
#include "EndomorphismSLP.h"

namespace crag {
namespace slp {
namespace recompression {
namespace {

std::string print_tree_preorder_single(const Vertex& vertex) {
  std::ostringstream out;
  std::unordered_set<slp::Vertex> printed;
  auto acceptor = [&printed] (const inspector::InspectorTask& task) {
    return printed.count(task.vertex) == 0;
    //true only if vertex is not visited yet and it is not terminal
  };
  Inspector<inspector::Preorder, decltype(acceptor)> inspector(vertex, acceptor);

  while (!inspector.stopped()) {
    PrintTo(inspector.vertex(), &out);
    out << " << ";
    if (printed.count(inspector.vertex().left_child()) || printed.count(inspector.vertex().left_child().negate())) {
      out << "(p) ";
    }
    PrintTo(inspector.vertex().left_child(), &out);
    out << " >> ";
    if (printed.count(inspector.vertex().right_child()) || printed.count(inspector.vertex().right_child().negate())) {
      out << "(p) ";
    }
    PrintTo(inspector.vertex().right_child(), &out);
    out << std::endl;
    printed.insert(inspector.vertex());

    ++inspector;
  }

  return out.str();
}

TEST(MadSorts, BitReverse) {
  ASSERT_EQ(0x8000000000000000ULL, mad_sorts::reverse_bits(0x0000000000000001ULL));
  ASSERT_EQ(0x0000000000000001ULL, mad_sorts::reverse_bits(0x8000000000000000ULL));
  ASSERT_EQ(0x0000000000000000ULL, mad_sorts::reverse_bits(0x0000000000000000ULL));
  ASSERT_EQ(0x0808080808080808ULL, mad_sorts::reverse_bits(0x1010101010101010ULL));
}

TEST(Recompression, ConstructionFromSLP) {
  TerminalVertex a(1);
  TerminalVertex b(2);
  NonterminalVertex ab(a, b);
  NonterminalVertex aa(a, a);
  JezRules rules(ab);
  JezRules rules_power(aa);

  ASSERT_EQ(1, rules.vertex_rules_.size());
  ASSERT_TRUE(rules.vertex_rules_.count(ab));

  Rule& rule = *(rules.vertex_rules_.begin()->second);

  EXPECT_EQ(2, rule.size());

  RuleLetter& a_node = *(rule.begin());
  RuleLetter& b_node = *(std::next(rule.begin()));
  EXPECT_TRUE(!a_node.is_nonterminal());
  EXPECT_TRUE(!a_node.is_power());
  EXPECT_TRUE(!b_node.is_nonterminal());
  EXPECT_TRUE(!b_node.is_power());
  EXPECT_EQ(1, a_node.terminal_power());
  EXPECT_EQ(1, b_node.terminal_power());
  EXPECT_TRUE(a_node.is_valid());
  EXPECT_TRUE(b_node.is_valid());

  EXPECT_GT(a_node.terminal_id(), 0);
  EXPECT_GT(b_node.terminal_id(), 0);
  EXPECT_EQ(3, a_node.terminal_id() + b_node.terminal_id());
  EXPECT_EQ(a_node.terminal_id(), rule.first_terminal_id());
  EXPECT_EQ(b_node.terminal_id(), rule.last_terminal_id());

  ASSERT_EQ(1, rules_power.vertex_rules_.size());
  ASSERT_TRUE(rules_power.vertex_rules_.count(aa));

  rule.delete_letter(rule.begin());
  EXPECT_TRUE(!a_node.is_valid());
  Rule::collect_garbage();
  rule.delete_letter(rule.begin());
  EXPECT_TRUE(rule.empty());

  Rule& rule_power = *(rules_power.vertex_rules_.begin()->second);

  EXPECT_EQ(1, rule_power.size());

  RuleLetter& node = *(rule_power.begin());
  EXPECT_TRUE(!node.is_nonterminal());
  EXPECT_TRUE(node.is_power());
  EXPECT_EQ(2, node.terminal_power());
  EXPECT_TRUE(node.is_valid());
  EXPECT_EQ(1, rule_power.first_terminal_id());
  EXPECT_EQ(1, rule_power.last_terminal_id());

  EXPECT_EQ(1, node.terminal_id());
}

TEST(Recompression, NormalForm) {
  TerminalVertex a(1);
  TerminalVertex b(2);
  NonterminalVertex ab(a, b);
  NonterminalVertex aa(a, a);
  NonterminalVertex aab(a, ab);
  NonterminalVertex aabaa(aab, aa);
  NonterminalVertex aba(ab, a);
  NonterminalVertex aaba(a, aba);

  Vertex normal_form_ab = normal_form(ab);
  EXPECT_EQ(a, normal_form_ab.left_child());
  EXPECT_EQ(b, normal_form_ab.right_child());
  EXPECT_EQ(VertexWord(ab), VertexWord(normal_form_ab));

  Vertex normal_form_aa = normal_form(aa);
  EXPECT_EQ(a, normal_form_aa.left_child());
  EXPECT_EQ(a, normal_form_aa.right_child());
  EXPECT_EQ(VertexWord(aa), VertexWord(normal_form_aa));

  Vertex normal_form_aab = normal_form(aab);
  EXPECT_EQ(a, normal_form_aab.left_child().left_child());
  EXPECT_EQ(a, normal_form_aab.left_child().right_child());
  EXPECT_EQ(b, normal_form_aab.right_child());
  EXPECT_EQ(VertexWord(aab), VertexWord(normal_form_aab));

  Vertex normal_form_aabaa = normal_form(aabaa);
  EXPECT_EQ(normal_form_aabaa.right_child().right_child(), normal_form_aabaa.left_child());
  EXPECT_EQ(a, normal_form_aabaa.left_child().left_child());
  EXPECT_EQ(a, normal_form_aabaa.left_child().right_child());
  EXPECT_EQ(b, normal_form_aabaa.right_child().left_child());
  EXPECT_EQ(VertexWord(aabaa), VertexWord(normal_form_aabaa));

  Vertex normal_form_aaba = normal_form(aaba);
  EXPECT_EQ(VertexWord(aaba), VertexWord(normal_form_aaba));
}

Vertex vertex_power(Vertex base, size_t power) {
  Vertex res;
  while (power) {
    if (power & 1) {
      if (res) {
        res = NonterminalVertex(res, base);
      } else {
        res = base;
      }
    }

    base = NonterminalVertex(base, base);
    power >>= 1;
  }

  return res;
}

namespace naive_Jez {

std::tuple<std::set<int>, std::set<int>>
greedy_pairs(std::set<std::pair<int, int>> pairs) {
  std::set<int> left_letters;
  std::set<int> right_letters;

  std::multimap<int, int> left_pairs;
  std::multimap<int, int> right_pairs;

  for (auto & pair : pairs) {
    left_pairs.insert(std::make_pair(pair.second, pair.first));
    right_pairs.insert(pair);
  }

  auto left_pair = left_pairs.begin();
  auto right_pair = right_pairs.begin();

  while (left_pair != left_pairs.end() || right_pair != right_pairs.end()) {
    int id = 0;
    if (left_pair == left_pairs.end()) {
      id = right_pair->first;
    } else if (right_pair == right_pairs.end()) {
      id = left_pair->first;
    } else {
      id = std::min(left_pair->first, right_pair->first);
    }

    decltype(left_pair) left_range_begin, left_range_end;
    decltype(right_pair) right_range_begin, right_range_end;

    size_t new_pairs_if_left = 0;

    if (right_pair->first == id) {

      std::tie(right_range_begin, right_range_end) =
          right_pairs.equal_range(right_pair->first);

      for (auto letter = right_range_begin; letter != right_range_end; ++letter) {
        if (letter->second > id || right_letters.count(letter->second)) {
          ++new_pairs_if_left;
        }
      }

      right_pair = right_range_end;
    }

    size_t new_pairs_if_right = 0;

    if (left_pair->first == id) {

      std::tie(left_range_begin, left_range_end) =
          left_pairs.equal_range(left_pair->first);

      for (auto letter = left_range_begin; letter != left_range_end; ++letter) {
        if (letter->second > id || left_letters.count(letter->second)) {
          ++new_pairs_if_right;
        }
      }
      left_pair = left_range_end;
    }

    if (new_pairs_if_left >= new_pairs_if_right) {
      left_letters.insert(id);
    } else {
      right_letters.insert(id);
    }
  }

  return std::make_tuple(std::move(left_letters), std::move(right_letters));
}

std::tuple<std::vector<Vertex>, std::vector<int>>
get_initial_rule(const Vertex& slp) {
  VertexWord vertex_word(slp);

  std::vector<Vertex> terminal_vertex;
  terminal_vertex.emplace_back();
  std::vector<int> word;

  std::map<int, int> letter_to_terminal;
  //it turns out that we reenumerate all terminals during normalization,
  //in order of their appearance in word.
  for (auto& letter : vertex_word) {
    if (!letter_to_terminal.count(letter)) {
      letter_to_terminal.insert(std::make_pair(letter, terminal_vertex.size()));
      terminal_vertex.push_back(TerminalVertex(letter));
    }
    word.push_back(letter_to_terminal[letter]);
  }

  return std::make_tuple(std::move(terminal_vertex), std::move(word));
}

std::set<std::pair<int, int>> get_pairs_list(const std::vector<int>& word) {
  std::set<std::pair<int, int>> pairs;
  auto current = word.begin();
  auto next = std::next(current);
  while (next != word.end()) {
    if (*current != *next) {
      pairs.insert(std::make_pair(*current, *next));
    }
    ++current;
    ++next;
  }

  return pairs;
}

bool compare_blocks(const std::pair<int, size_t>& first, const std::pair<int, size_t>& second) {
    return first.first < second.first ||
        (first.first == second.first && mad_sorts::reverse_bits(first.second) < mad_sorts::reverse_bits(second.second));
}

std::vector<int> compress_blocks(
    const std::vector<int>& word,
    std::vector<Vertex>* terminal_vertex
) {
  std::vector<std::pair<int, size_t>> blocks;

  for (auto& letter : word) {
    if (blocks.empty() || blocks.back().first != letter) {
      blocks.emplace_back(letter, 0);
    }
    ++blocks.back().second;
  }

  typedef bool(*PairsCompare)(const std::pair<int, size_t>&, const std::pair<int, size_t>&);
  std::map<std::pair<int, size_t>, int, PairsCompare> block_terminal(compare_blocks);

  for (auto& block : blocks) {
    if (block_terminal.count(block) == 0) {
      if (block.second > 1) {
        block_terminal[block] = 0;
      } else {
        block_terminal[block] = block.first;
      }
    }
  }

  for (auto& block : block_terminal) {
    if (!block.second) {
      Vertex base_terminal = terminal_vertex->at(block.first.first);
      block.second = terminal_vertex->size();
      terminal_vertex->push_back(vertex_power(
          base_terminal,
          block.first.second
      ));
    }
  }

  auto current_block = blocks.begin();
  auto current_block_terminal = block_terminal[*current_block];
  std::vector<int> without_blocks;
  for (auto& letter : word) {
    assert(current_block->first == letter);
    --current_block->second;
    if (current_block->second == 0) {
      without_blocks.push_back(current_block_terminal);
      ++current_block;
      current_block_terminal = block_terminal[*current_block];
    }
  }
  return without_blocks;
}

std::vector<int> compress_pairs(
    std::vector<int> word,
    const std::set<int>& left_letters,
    const std::set<int>& right_letters,
    std::set<std::pair<int, int>>* pairs,
    std::vector<Vertex>* terminal_vertex
) {
  if (word.size() < 2) {
    return word;
  }

  auto current_pair = pairs->begin();
  while (current_pair != pairs->end() && word.size() > 1) {
    auto next_pair = std::next(current_pair);

    if (left_letters.count(current_pair->first) &&
        right_letters.count(current_pair->second)) {

      NonterminalVertex pair_vertex(
          terminal_vertex->at(current_pair->first),
          terminal_vertex->at(current_pair->second)
      );

      std::vector<int> new_word;

      auto current = word.begin();
      auto next = std::next(current);

      while (next != word.end()) {
        if (*current == current_pair->first &&
            *next == current_pair->second) {
          new_word.push_back(terminal_vertex->size());
          current += 2;
          if (current == word.end()) {
            break;
          }
          next += 2;
        } else {
          new_word.push_back(*current);
          ++current;
          ++next;
        }
      }
      if (current != word.end()) {
        new_word.push_back(*current);
      }

      word = std::move(new_word);

      pairs->erase(current_pair);
      terminal_vertex->push_back(pair_vertex);
    }

    current_pair = next_pair;
  }

  return word;
}

Vertex normal_form(const Vertex& slp) {
  VertexWord vertex_word(slp);

  std::vector<Vertex> terminal_vertex;
  std::vector<int> word;

  std::tie(terminal_vertex, word) = get_initial_rule(slp);

  while (word.size() > 1) {
    std::set<std::pair<int, int>> pairs = get_pairs_list(word);
    word = compress_blocks(word, &terminal_vertex);

    while (!pairs.empty() && word.size() > 1) {
      std::set<int> left_letters;
      std::set<int> right_letters;

      std::tie(left_letters, right_letters) = greedy_pairs(pairs);

      word = compress_pairs(
          std::move(word),
          left_letters,
          right_letters,
          &pairs,
          &terminal_vertex
      );

      word = compress_pairs(
          std::move(word),
          right_letters,
          left_letters,
          &pairs,
          &terminal_vertex
      );

    }
  }

  return terminal_vertex[word.front()];
}

void debug_print_exposed(const std::vector<int>& word, ::std::ostream* os, int shift) {
  for (const int& letter : word) {
    (*os) << letter + shift << ',';
  }
}
} //namespace naive_Jez

::testing::AssertionResult is_normal_form(Vertex slp, Vertex normal_form) {
  if (VertexWord(slp) != VertexWord(normal_form)) {
    return ::testing::AssertionFailure() <<
        "Represented words are different" << std::endl <<
        VertexWord(slp) << " != " << VertexWord(normal_form) << std::endl;
  }
  Vertex naive_normal = naive_Jez::normal_form(slp);

  PostorderInspector first_inspector(normal_form);
  PostorderInspector second_inspector(naive_normal);

  while (!first_inspector.stopped()) {
    if (second_inspector.stopped()) {
      return testing::AssertionFailure();
    }

    VertexWord correct_word(second_inspector.vertex());
    VertexWord computed_word(first_inspector.vertex());

    if (correct_word != computed_word) {
      return ::testing::AssertionFailure() << "Different normal forms" <<
          "\n" << computed_word << " should be " << correct_word << "\n"
          << " for vertices "
          << ::testing::PrintToString(first_inspector.vertex())
          << " and "
          << ::testing::PrintToString(second_inspector.vertex())
          << "\n\nSlps:"
          << print_tree_preorder_single(normal_form) << "\n"
          << print_tree_preorder_single(naive_normal);
    }

    first_inspector.next();
    second_inspector.next();
  }

  return ::testing::AssertionSuccess();
}

::testing::AssertionResult normal_naive_equal(
    const Rule& rule,
    const std::vector<int>& word,
    int naive_shift
  ) {
  std::ostringstream normal_as_string;
  rule.debug_print_exposed(&normal_as_string);

  std::ostringstream naive_as_string;
  naive_Jez::debug_print_exposed(
      word,
      &naive_as_string,
      naive_shift
  );

  if (normal_as_string.str() == naive_as_string.str()) {
    return ::testing::AssertionSuccess();
  } else {
    return ::testing::AssertionFailure() <<
        normal_as_string.str() << " != " <<
        naive_as_string.str();
  }
}

//#define DEBUG_OUTPUT
void normalization_steps_check(const Vertex& root) {
  if (root.height() < 2) {
    return;
  }

  JezRules rules(root);

  Rule& root_rule = *(rules.vertex_rules_[root]);

  std::vector<Vertex> naive_terminal_vertex;
  std::vector<int> naive_word;

  std::tie(naive_terminal_vertex, naive_word) =
      naive_Jez::get_initial_rule(root);

  ASSERT_TRUE(normal_naive_equal(root_rule, naive_word, 0)) <<
      "Initial representations are not equal";

  while (root_rule.size() > 1 || (
           root_rule.size() == 1 && (
             root_rule.begin()->is_power() ||
             root_rule.begin()->is_nonterminal()
           )
        )) {
#ifdef DEBUG_OUTPUT
    std::cout << "\n=================\n\nCurrent rules:" << std::endl;
    rules.debug_print(&std::cout);
#endif
    OneStepPairs pairs(&rules);

    rules.remove_crossing_blocks();

#ifdef DEBUG_OUTPUT
    std::cout << "Rules after RemCrBlocks: " << std::endl;
    rules.debug_print(&std::cout);
#endif

    auto blocks = rules.list_blocks();
#ifdef DEBUG_OUTPUT
    std::cout << "\nFound blocks: " << std::endl;
    for (auto& block : blocks) {
      std::cout << block.rule_->debug_id << ':';
      block.letter_->debug_print(&std::cout);
      std::cout << std::endl;
    }
#endif

    rules.compress_blocks(blocks);

#ifdef DEBUG_OUTPUT
    std::cout << "Rules after CompressBlocks: " << std::endl;
    rules.debug_print(&std::cout);
#endif

    ASSERT_GE(naive_word.size(), 2);

    std::set<std::pair<int, int>> naive_pairs =
        naive_Jez::get_pairs_list(naive_word);

    naive_word = naive_Jez::compress_blocks(
        naive_word,
        &naive_terminal_vertex
    );

    ASSERT_TRUE(normal_naive_equal(root_rule, naive_word, 0)) <<
        "Representations after blocks compression are not equal";

    std::vector<unsigned char> left_letters, right_letters;

    std::tie(left_letters, right_letters) = pairs.greedy_pairs();
#ifdef DEBUG_OUTPUT
    std::cout << "\nGreedyPairs:\nLeft: ";

    for (auto& terminal : left_letters) {
      std::cout << terminal << ',';
    }
    std::cout << "\nRight: ";
    for (auto& terminal : right_letters) {
      std::cout << terminal << ',';
    }
    std::cout << std::endl;
#endif

    std::set<int> naive_left_letters;
    std::set<int> naive_right_letters;

    while (!left_letters.empty()) {
      ASSERT_TRUE(!naive_pairs.empty());

      std::tie(naive_left_letters, naive_right_letters) =
          naive_Jez::greedy_pairs(naive_pairs);

      ASSERT_EQ(naive_left_letters.size(), count(left_letters.begin(), left_letters.end(), 1));
      ASSERT_EQ(naive_right_letters.size(), count(right_letters.begin(), right_letters.end(), 1));

      for (auto& letter : naive_left_letters) {
        ASSERT_TRUE(letter < left_letters.size() && left_letters.at(letter))
            << "Left letter " << letter <<
            " is not among normal left letters";
      }

      for (auto& letter : naive_right_letters) {
        ASSERT_TRUE(letter < right_letters.size() && right_letters.at(letter))
            << "Right letter " << letter
            << " is not among normal right letters";
      }

      pairs.remove_crossing(left_letters, right_letters);
      pairs.compress_pairs_from_letter_lists(left_letters, right_letters);
#ifdef DEBUG_OUTPUT
      std::cout << "Rules after first compression: " << std::endl;
      rules.debug_print(&std::cout);
#endif

      naive_word = naive_Jez::compress_pairs(
          std::move(naive_word),
          naive_left_letters,
          naive_right_letters,
          &naive_pairs,
          &naive_terminal_vertex
      );

      ASSERT_TRUE(normal_naive_equal(root_rule, naive_word, 0)) <<
          "Representations after first pairs compression are not equal";

      pairs.remove_crossing(right_letters, left_letters);
      pairs.compress_pairs_from_letter_lists(right_letters, left_letters);
#ifdef DEBUG_OUTPUT
      std::cout << "Rules after second compression: " << std::endl;
      rules.debug_print(&std::cout);
#endif

      naive_word = naive_Jez::compress_pairs(
          std::move(naive_word),
          naive_right_letters,
          naive_left_letters,
          &naive_pairs,
          &naive_terminal_vertex
      );

      ASSERT_TRUE(normal_naive_equal(root_rule, naive_word, 0)) <<
          "Representations after second pairs compression are not equal";


      std::tie(left_letters, right_letters) = pairs.greedy_pairs();
#ifdef DEBUG_OUTPUT
      std::cout << "\nGreedyPairs:\nLeft: ";

      for (auto& terminal : left_letters) {
        std::cout << terminal << ',';
      }
      std::cout << "\nRight: ";
      for (auto& terminal : right_letters) {
        std::cout << terminal << ',';
      }
      std::cout << std::endl;
#endif
    }

    ASSERT_TRUE(naive_pairs.empty());

    rules.empty_cleanup();
    Rule::collect_garbage();
  }
}



TEST(Recompression, NormalFormEx1) {
  TerminalVertex a(1);
  TerminalVertex b(2);
  NonterminalVertex ab(a, b);
  NonterminalVertex aba(ab, a);
  NonterminalVertex aaba(a, aba);

  EXPECT_TRUE(is_normal_form(aaba, normal_form(aaba)));
  normalization_steps_check(aaba);
}

TEST(Recompression, NormalFormEx2) {
  TerminalVertex a(1);
  TerminalVertex b(2);
  NonterminalVertex ab(a, b);
  NonterminalVertex ba(b, a);
  NonterminalVertex abba(ab, ba);
  NonterminalVertex babba(b, abba);

  EXPECT_TRUE(is_normal_form(babba, normal_form(babba)));
  normalization_steps_check(babba);
}


TEST(Recompression, NormalFormEx3) {
  TerminalVertex a(1);
  TerminalVertex b(2);
  TerminalVertex c(3);
  NonterminalVertex bc(b, c);
  NonterminalVertex bcb(bc, b);
  NonterminalVertex bbcb(b, bcb);
  NonterminalVertex bbcba(bbcb, a);
  NonterminalVertex abbcba(a, bbcba);

  Vertex normal_form_abbcba = normal_form(abbcba);
  EXPECT_TRUE(is_normal_form(abbcba, normal_form(abbcba)));
  normalization_steps_check(abbcba);
}

TEST(Recompression, NormalFormEx4) {
  TerminalVertex a(1);
  TerminalVertex b(2);
  NonterminalVertex ab(a, b);
  NonterminalVertex bab(b, ab);
  NonterminalVertex babab(bab, ab);

  EXPECT_TRUE(is_normal_form(babab, normal_form(babab)));
  normalization_steps_check(babab);
}

TEST(Recompression, NormalFormEx5) {
  TerminalVertex a(1);
  TerminalVertex b(2);
  NonterminalVertex ba(b, a);
  NonterminalVertex baa(ba, a);
  NonterminalVertex babaa(ba, baa);

  EXPECT_TRUE(is_normal_form(babaa, normal_form(babaa)));
  normalization_steps_check(babaa);
}

TEST(Recompression, NormalFormEx6) {
  TerminalVertex a(1);
  NonterminalVertex aa(a, a);
  NonterminalVertex aaa(aa, a);
  NonterminalVertex aaaa(a, aaa);

  EXPECT_TRUE(is_normal_form(aaaa, normal_form(aaaa)));
  normalization_steps_check(aaaa);
}

TEST(Recompression, NormalFormEx7) {
  TerminalVertex a(1);
  TerminalVertex b(2);
  TerminalVertex c(3);
  NonterminalVertex ba(b, a);
  NonterminalVertex aba(a, ba);
  NonterminalVertex baba(b, aba);
  NonterminalVertex cbaba(c, baba);

  EXPECT_TRUE(is_normal_form(cbaba, normal_form(cbaba)));
  normalization_steps_check(cbaba);
}


Vertex get_random_slp_on_n_letters(unsigned int WORD_SIZE, unsigned int ALPHABET_SIZE) {
  std::vector<unsigned int> random_word_split;
  for (unsigned int i = 1; i < WORD_SIZE; ++i) {
    random_word_split.push_back(i);
  }

  std::random_shuffle(random_word_split.begin(), random_word_split.end());

  std::vector<Vertex> word_presentation;
  for (unsigned int i = 0; i < WORD_SIZE; ++i) {
    word_presentation.push_back(TerminalVertex(rand() % ALPHABET_SIZE + 1));
  }

  for (unsigned int split : random_word_split) {
    NonterminalVertex new_vertex(word_presentation[split - 1], word_presentation[split]);
    for (unsigned int i = split - new_vertex.left_child().length().get_ui(); i < split + new_vertex.right_child().length(); ++i) {
      word_presentation[i] = new_vertex;
    }
  }

  return word_presentation.front();
}


TEST(Recompression, NormalFormEx8) {
  TerminalVertex a(1);
  TerminalVertex b(2);
  NonterminalVertex ab(a, b);
  NonterminalVertex bab(b, ab);
  NonterminalVertex babab(bab, ab);
  NonterminalVertex babbabab(bab, babab);

  EXPECT_TRUE(is_normal_form(babbabab, normal_form(babbabab)));
  normalization_steps_check(babbabab);
}

TEST(Recompression, NormalFormEx9) {
  TerminalVertex a(1);
  TerminalVertex b(2);
  TerminalVertex c(3);
  NonterminalVertex bc(b, c);
  NonterminalVertex bcb(bc, b);
  NonterminalVertex abcb(a, bcb);
  NonterminalVertex cabcb(c, abcb);

  EXPECT_TRUE(is_normal_form(cabcb, normal_form(cabcb)));
  normalization_steps_check(cabcb);
}

TEST(Recompression, NormalFormEx10) {
  TerminalVertex a(1);
  TerminalVertex b(2);
  NonterminalVertex ab(a, b);
  NonterminalVertex bab(b, ab);
  NonterminalVertex babab(bab, ab);
  NonterminalVertex bababbab(babab, bab);

  EXPECT_TRUE(is_normal_form(bababbab, normal_form(bababbab)));
  normalization_steps_check(bababbab);
}

TEST(Recompression, NormalFormEx11) {
  TerminalVertex a(1);
  TerminalVertex b(2);
  NonterminalVertex ab(a, b);
  NonterminalVertex aab(a, ab);
  NonterminalVertex aaab(a, aab);
  NonterminalVertex aaabaab(aaab, aab);

  //EXPECT_TRUE(is_normal_form(aaabaab, normal_form(aaabaab)));
  normalization_steps_check(aaabaab);
}

TEST(Recompression, NormalFormEx12) {
  TerminalVertex a(1);
  TerminalVertex b(2);
  NonterminalVertex ab(a, b);
  NonterminalVertex aab(a, ab);
  NonterminalVertex aabab(aab, ab);
  NonterminalVertex aababaab(aabab, aab);
  NonterminalVertex aababaabaabab(aababaab, aabab);
  NonterminalVertex aababaabaababaababaab(aababaabaabab, aababaab);

  //EXPECT_TRUE(is_normal_form(aababaabaababaababaab, normal_form(aababaabaababaababaab)));
  normalization_steps_check(aababaabaababaababaab);
}

TEST(Recompression, NormalFormEx13) {
  TerminalVertex a(1);
  TerminalVertex b(2);
  TerminalVertex c(3);
  NonterminalVertex ab(a, b);       //0
  NonterminalVertex cab(c, ab);     //1
  NonterminalVertex abcab(ab, cab); //2
  NonterminalVertex cb(c, b);       //3
  NonterminalVertex abcb(ab, cb);   //4
  NonterminalVertex abcababcb(abcab, abcb); //5
  NonterminalVertex abcababcababcb(abcab, abcababcb); //6

  EXPECT_TRUE(is_normal_form(abcababcababcb, normal_form(abcababcababcb)));
  normalization_steps_check(abcababcababcb);
}

TEST(Recompression, NormalFormEx14) {
  TerminalVertex a(1);
  TerminalVertex b(2);
  TerminalVertex c(3);
  NonterminalVertex v0(a, b); //ab
  NonterminalVertex v1(v0, a); //aba
  NonterminalVertex v2(c, v1); //caba
  NonterminalVertex v3(v0, v2); //abcaba
  NonterminalVertex v4(v3, v1); //abcabaaba
  NonterminalVertex v5(v1, v0); //abaab
  NonterminalVertex v6(v4, v5); //abcabaabaabaab
  NonterminalVertex v7(v4, v6); //abcabaabaabcabaabaabaab
  NonterminalVertex v8(v4, v7); //abcabaabaabcabaabaabcabaabaabaab
  NonterminalVertex v9(v1, v8); //abaabcabaabaabcabaabaabcabaabaabaab
  NonterminalVertex v10(v1, v4); //abaabcabaaba
  NonterminalVertex v11(v10, v4); //abaabcabaabaabcabaaba
  NonterminalVertex v12(v9, v11); //abaabcabaabaabcabaabaabcabaabaabaababaabcabaabaabcabaaba
  NonterminalVertex v13(v11, v4); //abaabcabaabaabcabaabaabcabaaba
  NonterminalVertex v14(v13, v4); //abaabcabaabaabcabaabaabcabaabaabcabaaba
  NonterminalVertex v15(v12, v14); //abaabcabaabaabcabaabaabcabaabaabaababaabcabaabaabcabaabaabaabcabaabaabcabaabaabcabaabaabcabaaba

  EXPECT_TRUE(is_normal_form(v15, normal_form(v15)));
  normalization_steps_check(v15);
}

TEST(Recompression, NormalFormEx15) {
  TerminalVertex a(1);
  TerminalVertex b(2);
  TerminalVertex c(3);
  NonterminalVertex v0(a, b); //ab
  NonterminalVertex v1(v0, c); //abc
  NonterminalVertex v2(v1, v0); //abcab
  NonterminalVertex v3(v2, v0); //abcabab
  NonterminalVertex v4(v0, v3); //ababcabab
  NonterminalVertex v5(v4, v0); //ababcababab
  NonterminalVertex v6(v5, v0); //ababcabababab
  NonterminalVertex v7(c, b); //cb
  NonterminalVertex v8(v7, v5); //cbababcababab
  NonterminalVertex v9(v0, v8); //abcbababcababab
  NonterminalVertex v10(v6, v9); //ababcabababababcbababcababab
  NonterminalVertex v11(v0, v6); //abababcabababab
  NonterminalVertex v12(v0, v11); //ababababcabababab
  NonterminalVertex v13(v0, v12); //abababababcabababab
  NonterminalVertex v14(v13, v12); //abababababcababababababababcabababab
  NonterminalVertex v15(v14, v12); //abababababcababababababababcababababababababcabababab
  NonterminalVertex v16(v10, v15); //ababcabababababcbababcababababababababcababababababababcababababababababcabababab

  EXPECT_TRUE(is_normal_form(v16, normal_form(v16)));
  normalization_steps_check(v16);
}

TEST(Recompression, NormalFormEx16) {
  TerminalVertex a(1);
  TerminalVertex b(2);
  TerminalVertex c(3);
  NonterminalVertex v0(a, b); //ab
  NonterminalVertex v1(b, c); //bc
  NonterminalVertex v2(a, v1); //abc
  NonterminalVertex v3(v1, v2); //bcabc
  NonterminalVertex v4(v0, v3); //abbcabc
  NonterminalVertex v5(v4, v2); //abbcabcabc
  NonterminalVertex v6(v5, v2); //abbcabcabcabc
  NonterminalVertex v7(v6, v5); //abbcabcabcabcabbcabcabc
  NonterminalVertex v8(v7, v6); //abbcabcabcabcabbcabcabcabbcabcabcabc
  NonterminalVertex v9(v6, v3); //abbcabcabcabcbcabc
  NonterminalVertex v10(v6, v9); //abbcabcabcabcabbcabcabcabcbcabc
  NonterminalVertex v11(v8, v10); //abbcabcabcabcabbcabcabcabbcabcabcabcabbcabcabcabcabbcabcabcabcbcabc
  NonterminalVertex v12(v6, v10); //abbcabcabcabcabbcabcabcabcabbcabcabcabcbcabc
  NonterminalVertex v13(v11, v12); //abbcabcabcabcabbcabcabcabbcabcabcabcabbcabcabcabcabbcabcabcabcbcabcabbcabcabcabcabbcabcabcabcabbcabcabcabcbcabc
  NonterminalVertex v14(v13, v10); //abbcabcabcabcabbcabcabcabbcabcabcabcabbcabcabcabcabbcabcabcabcbcabcabbcabcabcabcabbcabcabcabcabbcabcabcabcbcabcabbcabcabcabcabbcabcabcabcbcabc
  NonterminalVertex v15(v14, v13); //abbcabcabcabcabbcabcabcabbcabcabcabcabbcabcabcabcabbcabcabcabcbcabcabbcabcabcabcabbcabcabcabcabbcabcabcabcbcabcabbcabcabcabcabbcabcabcabcbcabcabbcabcabcabcabbcabcabcabbcabcabcabcabbcabcabcabcabbcabcabcabcbcabcabbcabcabcabcabbcabcabcabcabbcabcabcabcbcabc

  normalization_steps_check(v15);
  EXPECT_TRUE(is_normal_form(v15, normal_form(v15)));
}

TEST(Recompression, NormalFormEx17) {
  TerminalVertex a(1);
  TerminalVertex b(2);
  TerminalVertex c(3);
  NonterminalVertex v0(a, b); //
  NonterminalVertex v1(v0, a); //
  NonterminalVertex v2(v0, v1); //
  NonterminalVertex v3(v1, v2); //
  NonterminalVertex v4(v2, v3); //
  NonterminalVertex v5(v2, v4); //
  NonterminalVertex v6(v4, v5); //
  NonterminalVertex v7(v5, v6); //
  NonterminalVertex v8(v6, v7); //
  NonterminalVertex v9(v8, v7); //
  NonterminalVertex v10(v5, c); //
  NonterminalVertex v11(v6, v10); //
  NonterminalVertex v12(v6, v11); //
  NonterminalVertex v13(v9, v12); //
  NonterminalVertex v14(v7, v13); //
  NonterminalVertex v15(v13, v14); //
  NonterminalVertex v16(v13, v15); //
  NonterminalVertex v17(v13, v16); //
  NonterminalVertex v18(v16, v17); //

  normalization_steps_check(v18);
}

TEST(Recompression, NormalFormEx18) {
  TerminalVertex a(1);
  TerminalVertex b(2);
  TerminalVertex c(3);
  NonterminalVertex v0(a, c); //
  NonterminalVertex v1(b, v0); //
  NonterminalVertex v2(a, v1); //
  NonterminalVertex v3(v1, v0); //
  NonterminalVertex v4(v2, v3); //

  normalization_steps_check(v4);
}

TEST(Recompression, NormalFormEx19) {
  TerminalVertex a(1);
  TerminalVertex b(2);
  TerminalVertex c(3);
  NonterminalVertex v0(a, b); //
  NonterminalVertex v1(c, b); //
  NonterminalVertex v2(v0, v1); //
  NonterminalVertex v3(v2, v1); //

  normalization_steps_check(v3);
}


std::string print_rules(const Vertex& slp) {
  char fresh_terminal = 'a' - 1;
  std::unordered_map<Vertex, std::string> vertex_rules;
  //std::unordered_map<Vertex, std::string> vertex_strings;
  std::stringstream rules;
  auto acceptor = [&vertex_rules] (const inspector::InspectorTask& task) {
    return vertex_rules.count(task.vertex) == 0;
    //true only if vertex is not visited yet and it is not terminal
  };
  Inspector<inspector::Postorder, decltype(acceptor)> inspector(slp, acceptor);
  int rule_count = 0;
  while (!inspector.stopped()) {
    if (inspector.vertex().height() < 2) {
      auto inserted = vertex_rules.emplace(
          inspector.vertex(),
          std::string(1, ++fresh_terminal)
      );
//      vertex_strings[inspector.vertex()] = inserted.first->second;
    } else {
      Vertex left = inspector.vertex().left_child();
      Vertex right = inspector.vertex().right_child();

      std::string& left_rule = vertex_rules[left];
      std::string& right_rule = vertex_rules[right];
      std::stringstream vertex_rule;
      vertex_rule << "v" << rule_count++;
      vertex_rules[inspector.vertex()] = vertex_rule.str();
//      vertex_strings[inspector.vertex()] = vertex_strings[left] + vertex_strings[right];
      rules << "  NonterminalVertex "
            << vertex_rules[inspector.vertex()]
            << "(" << left_rule << ", " << right_rule << "); //"
//            << vertex_strings[inspector.vertex()]
            << "\n";
    }

    inspector.next();
  }
  std::string terminal_rules;
  for (char c = 'a'; c <= fresh_terminal; ++c) {
    terminal_rules += "  TerminalVertex ";
    terminal_rules += c;
    terminal_rules += "(";
    terminal_rules += c - 'a' + '1';
    terminal_rules += ");\n";
  }

  rules << "\n  Vertex normal_form_";
  rules << vertex_rules[slp];
  rules << "= normal_form(";
  rules << vertex_rules[slp];
  rules << ");\n";

  rules << "EXPECT_TRUE(is_normal_form(";
  rules << vertex_rules[slp];
  rules << ", normal_form(";
  rules << vertex_rules[slp];
  rules << ")));\n";

  return terminal_rules + rules.str();
}

TEST(Recompression, StressNormalForm) {
  const unsigned int WORD_SIZE = 100;
  const unsigned int ALPHABET_SIZE = 3;
  int REPEAT = 1000;

  typedef crag::slp::TVertexHashAlgorithms<
      crag::slp::hashers::SinglePowerHash,
      crag::slp::hashers::PermutationHash<crag::Permutation16>
  > VertexHashAlgorithms;

  srand(1717);

  while (--REPEAT >= 0) {
    Vertex slp = get_random_slp_on_n_letters(WORD_SIZE, ALPHABET_SIZE);
    //std::cout << print_rules(slp) << std::endl;
    Vertex normal_slp = normal_form(slp);

    ASSERT_EQ(
        slp.length(),
        normal_slp.length()
    );

    VertexHashAlgorithms::Cache hash_cache;

    ASSERT_EQ(
        VertexHashAlgorithms::get_subvertex_hash(slp, 0, slp.length()),
        VertexHashAlgorithms::get_subvertex_hash(normal_slp, 0, normal_slp.length(), &hash_cache)
    ) << print_tree_preorder_single(slp) << "\n"
      << print_tree_preorder_single(normal_slp);


    std::unordered_set<slp::Vertex> checked;
    std::unordered_set<VertexHashAlgorithms::VertexHash> hashes;

    auto acceptor = [&checked] (const inspector::InspectorTask& task) {
      return checked.count(task.vertex) == 0;
      //true only if vertex is not visited yet and it is not terminal
    };
    Inspector<inspector::Preorder, decltype(acceptor)> inspector(normal_slp, acceptor);

    while (!inspector.stopped()) {
      auto hash = VertexHashAlgorithms::get_subvertex_hash(
          inspector.vertex(),
          0,
          inspector.vertex().length(),
          &hash_cache
      );

      auto inserted = hashes.insert(hash);

      ASSERT_TRUE(inserted.second) << VertexWord(inspector.vertex());

//      if (!inserted.second) {
//        VertexWord word(inspector.vertex());
//        auto terminal = word[0];
//
//        for (auto& symbol : word) {
//          ASSERT_EQ(terminal, symbol);
//        }
//      }
//
//      ASSERT_TRUE(inserted.second) << print_tree_preorder_single(slp) << "\n"
//          << print_tree_preorder_single(normal_slp);


      checked.insert(inspector.vertex());
      inspector.next();
    }

  }
}

TEST(Recompression, StressEndomorphismNormal) {
  CONSTEXPR_OR_CONST size_t RANK = 3;
  CONSTEXPR_OR_CONST size_t ENDOMORPHISMS_NUMBER = 50;
  int REPEAT = 1000;
  size_t seed = 112233;
  srand(seed);
  UniformAutomorphismSLPGenerator<> generator(RANK, seed);
  generator.set_inverters_probability(0);

  typedef crag::slp::TVertexHashAlgorithms<
      crag::slp::hashers::SinglePowerHash,
      crag::slp::hashers::PermutationHash<crag::Permutation16>
  > VertexHashAlgorithms;

  while (--REPEAT >= 0) {
    Vertex slp = EndomorphismSLP::composition(ENDOMORPHISMS_NUMBER, generator).image(1);
    //std::cout << print_rules(slp) << std::endl;
    Vertex normal_slp = normal_form(slp);

    ASSERT_EQ(
        slp.length(),
        normal_slp.length()
    );

    VertexHashAlgorithms::Cache hash_cache;

    ASSERT_EQ(
        VertexHashAlgorithms::get_subvertex_hash(slp, 0, slp.length()),
        VertexHashAlgorithms::get_subvertex_hash(normal_slp, 0, normal_slp.length(), &hash_cache)
    ) << print_tree_preorder_single(slp) << "\n"
      << print_tree_preorder_single(normal_slp) << "\n"
      << print_rules(slp);


    std::unordered_set<slp::Vertex> checked;
    std::unordered_map<VertexHashAlgorithms::VertexHash, Vertex> hashes;

    auto acceptor = [&checked] (const inspector::InspectorTask& task) {
      return checked.count(task.vertex) == 0;
      //true only if vertex is not visited yet and it is not terminal
    };
    Inspector<inspector::Preorder, decltype(acceptor)> inspector(normal_slp, acceptor);

    while (!inspector.stopped()) {
      auto hash = VertexHashAlgorithms::get_subvertex_hash(
          inspector.vertex(),
          0,
          inspector.vertex().length(),
          &hash_cache
      );

      auto inserted = hashes.emplace(hash, inspector.vertex());

      ASSERT_TRUE(inserted.second) << print_tree_preorder_single(slp) << "\n"
          << print_tree_preorder_single(normal_slp) << "\n"
          << print_rules(slp) << "\n\n"
          << print_tree_preorder_single(inserted.first->second) << "\n"
          << print_tree_preorder_single(inspector.vertex()) << "\n"
          << VertexWord(inspector.vertex());
//      if (!inserted.second) {
//        VertexWord word(inspector.vertex());
//        auto terminal = word[0];
//
//        for (auto& symbol : word) {
//          ASSERT_EQ(terminal, symbol) << print_tree_preorder_single(slp) << "\n"
//              << print_tree_preorder_single(normal_slp) << "\n"
//              << print_rules(slp) << "\n\n"
//              << print_tree_preorder_single(inserted.first->second) << "\n"
//              << print_tree_preorder_single(inspector.vertex()) << "\n"
//              << VertexWord(inspector.vertex());
//        }
//      }

//      ASSERT_TRUE(inserted.second) << print_tree_preorder_single(slp) << "\n"
//          << print_tree_preorder_single(normal_slp);


      checked.insert(inspector.vertex());
      inspector.next();
    }

  }

}



} //namespace
}
} //slp
} //crag


