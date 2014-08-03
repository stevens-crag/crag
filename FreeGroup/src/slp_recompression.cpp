#include "slp_recompression.h"
#include <iterator>

namespace crag {
namespace slp {
namespace recompression {

Rule::Rule(std::initializer_list<RuleLetter> letters)
  : debug_id(0)
{
  for (auto& letter : letters) {
    if (!letters_.empty() &&
        !letter.is_nonterminal() &&
        letters_.back().is_power_of(letter.terminal_id())) {
      letters_.back().terminal_.power += letter.terminal_power();
    } else {
      letters_.push_back(letter);

      if (letter.is_nonterminal()) {
        letters_.back().nonterminal_rule_->register_inclusion(
            this,
            std::prev(letters_.end())
        );
      }
    }
  }

  if (!letters_.empty()) {
    first_terminal_letter_ = letters_.front().first_terminal_ptr();
    last_terminal_letter_ = letters_.back().last_terminal_ptr();
  }

  assert(first_terminal_letter_ && last_terminal_letter_);
}

Rule::iterator Rule::pop_first_from_letter(Rule::iterator letter_position) {
  assert(!empty());
  assert(!letter_position->is_empty_nonterminal());

  if (!letter_position->is_nonterminal()) {
    return letter_position;
  }

  Rule* letter_rule = letter_position->nonterminal_rule();

  iterator popping_letter = letter_rule->begin();
  assert(!popping_letter->is_empty_nonterminal());

  if (popping_letter->is_nonterminal()) {
    popping_letter = letter_rule->pop_first_from_letter(popping_letter);
  }

  assert(!popping_letter->is_nonterminal());
  assert(popping_letter == letter_rule->begin());

  const RuleLetter& popped_letter = *(letter_rule->delete_letter(popping_letter));

  assert(letter_rule->empty() || !letter_rule->begin()->is_empty_nonterminal());

  if (!letter_rule->empty()) {
    letter_rule->first_terminal_letter_ = letter_rule->begin()->first_terminal_ptr();
  } else {
    letter_rule->first_terminal_letter_ = nullptr;
    letter_rule->first_terminal_letter_ = nullptr;
  }

  iterator inserted_letter;

  for (auto& occurence : letter_rule->nonterminal_index_) {
    assert(!occurence.rule_->empty());

    occurence.rule_->insert_popped_letter_left(occurence.letter_, popped_letter);

    if (occurence.rule_ == this && occurence.letter_ == letter_position) {
      if (letter_rule->empty()) {
        assert(occurence.letter_->is_empty_nonterminal());
        inserted_letter = occurence.rule_->remove_empty_letter(occurence.letter_).first;
        assert(inserted_letter != end());
      } else {
        assert(letter_position != begin());
        inserted_letter = std::prev(letter_position);
      }
    } else {
      if (letter_rule->empty()) {
        assert(occurence.letter_->is_empty_nonterminal());
        occurence.rule_->remove_empty_letter(occurence.letter_);
      }
    }
  }

  return inserted_letter;
}

Rule::iterator Rule::pop_last_from_letter(Rule::iterator letter_position) {
  assert(!empty());
  assert(!letter_position->is_empty_nonterminal());

  if (!letter_position->is_nonterminal()) {
    return letter_position;
  }

  Rule* letter_rule = letter_position->nonterminal_rule();

  iterator popping_letter = std::prev(letter_rule->end());
  assert(!popping_letter->is_empty_nonterminal());

  if (popping_letter->is_nonterminal()) {
    popping_letter = letter_rule->pop_last_from_letter(popping_letter);
  }

  assert(!popping_letter->is_nonterminal());
  assert(popping_letter == std::prev(letter_rule->end()));

  const RuleLetter& popped_letter = *(letter_rule->delete_letter(popping_letter));

  assert(letter_rule->empty() || !std::prev(letter_rule->end())->is_empty_nonterminal());

  if (!letter_rule->empty()) {
    letter_rule->last_terminal_letter_ = std::prev(letter_rule->end())->last_terminal_ptr();
  } else {
    letter_rule->first_terminal_letter_ = nullptr;
    letter_rule->first_terminal_letter_ = nullptr;
  }

  iterator inserted_letter;

  for (auto& occurence : letter_rule->nonterminal_index_) {
    assert(!occurence.rule_->empty());

    occurence.rule_->insert_popped_letter_right(occurence.letter_, popped_letter);

    if (occurence.rule_ == this && occurence.letter_ == letter_position) {
      if (letter_rule->empty()) {
        assert(occurence.letter_->is_empty_nonterminal());
        inserted_letter = occurence.rule_->remove_empty_letter(occurence.letter_).second;
        assert(inserted_letter != end());
      } else {
        inserted_letter = std::next(letter_position);
        assert(inserted_letter != end());
      }
    } else {
      if (letter_rule->empty()) {
        assert(occurence.letter_->is_empty_nonterminal());
        occurence.rule_->remove_empty_letter(occurence.letter_);
      }
    }
  }

  return inserted_letter;
}

std::pair<Rule::iterator, Rule::iterator> Rule::remove_empty_letter(Rule::iterator position) {
  assert(!empty());
  assert(position->is_empty_nonterminal());

  auto position_after = std::next(position);
  auto position_before = std::prev(position);

  delete_letter(position);

  if (empty()) {
    first_terminal_letter_ = nullptr;
    last_terminal_letter_ = nullptr;
    return std::make_pair(end(), end());
  }

  if (position_after != letters_.begin() &&
      position_after != letters_.end() &&
      !position_before->is_nonterminal() &&
      position_after->is_power_of(position_before->terminal_id())) {
    position_before->terminal_.power += position_after->terminal_power();
    delete_letter(position_after);

    assert(!std::prev(end())->is_empty_nonterminal());

    return std::make_pair(position_before, position_before);
  }

  assert(!empty());

  if (position_after == letters_.begin()) {
    return std::make_pair(end(), std::move(position_after));
  } else if (position_before == std::prev(end())) {
    return std::make_pair(position_before, end());
  }

  return std::make_pair(std::move(position_before), std::move(position_after));
}

Rule::iterator Rule::compress_power(Rule::iterator position, TerminalId new_terminal) {
  assert(!position->is_nonterminal());

  assert(!begin()->is_empty_nonterminal());
  assert(!std::prev(end())->is_empty_nonterminal());

  position->terminal_.id = new_terminal;
  position->terminal_.power = 1;

  assert(position == begin() || !std::prev(position)->is_power_of(new_terminal));
  assert(std::next(position) == end() || !std::next(position)->is_power_of(new_terminal));

  return position;
}

Rule::iterator Rule::compress_pair(
    Rule::iterator first,
    Rule::iterator second,
    TerminalId new_terminal
) {
  assert(second == std::next(first));
  assert(!first->is_nonterminal() && !first->is_power());
  assert(!second->is_nonterminal() && !second->is_power());

  assert(!begin()->is_empty_nonterminal());

  assert(!std::prev(end())->is_empty_nonterminal());

  first->terminal_.id = new_terminal;
  delete_letter(second); //we have to delete the second letter after we have introduced new terminal

  auto prev = std::prev(first);

  if (first != letters_.begin() &&
      prev->is_power_of(new_terminal)) {
    prev->terminal_.power += 1;
    delete_letter(first);
    first = prev;
  }

  assert(std::next(first) == end() || !std::next(first)->is_power_of(new_terminal));

  return first;
}

void Rule::insert_popped_letter_right(
    Rule::iterator letter_position,
    const RuleLetter& popped_letter)
{
  assert(!popped_letter.is_nonterminal());
  assert(!popped_letter.is_valid());
  assert(letter_position->is_nonterminal());

  auto position_after = std::next(letter_position);

  if (position_after != end() &&
      position_after->is_power_of(popped_letter.terminal_id())) {

    position_after->terminal_.power += popped_letter.terminal_power();
    return;
  }

  letters_.emplace(position_after, popped_letter);

  //assert(!std::prev(end())->is_empty_nonterminal());
}

void Rule::insert_popped_letter_left(
    Rule::iterator letter_position,
    const RuleLetter& popped_letter)
{
  assert(!popped_letter.is_nonterminal());
  assert(!popped_letter.is_valid());
  assert(letter_position->is_nonterminal());

  auto position_before = std::prev(letter_position);
  if (letter_position != letters_.begin() &&
      position_before->is_power_of(popped_letter.terminal_id())) {

    position_before->terminal_.power += popped_letter.terminal_power();
    return;
  }

  letters_.emplace(letter_position, popped_letter);

  //assert(!begin()->is_empty_nonterminal());
}

JezRules::JezRules(const Vertex& slp)
  : fresh_terminal_id_(0)
{
  auto acceptor = [this] (const inspector::InspectorTask& task) {
    return this->vertex_rules_.count(task.vertex) == 0;
    //true only if vertex is not visited yet
  };
  Inspector<inspector::Postorder, decltype(acceptor)> inspector(slp, acceptor);

  while (!inspector.stopped()) {
    if (inspector.vertex().height() < 2) {
      get_letter(inspector.vertex());
    } else {
      Vertex left = inspector.vertex().left_child();
      Vertex right = inspector.vertex().right_child();

      rules_.emplace_back<std::initializer_list<RuleLetter>>({
        get_letter(left),
        get_letter(right)
      });
      rules_.back().debug_id = rules_.size() - 1;

      vertex_rules_.insert(
        std::make_pair(
          inspector.vertex(),
          &(rules_.back())
        )
      );
    }
    inspector.next();
  }
}

void JezRules::remove_crossing_blocks() {
  for (auto& rule : rules_) {
    auto current = rule.begin();


    while (current != rule.end()) {
      assert(!current->is_empty_nonterminal());
      assert(current->is_valid());
      assert(!current->is_empty_nonterminal());
      auto next = std::next(current);
      assert(next == rule.end() || !next->is_empty_nonterminal());

      if (next != rule.end() &&
          current->last_terminal_letter_id() ==
          next->first_terminal_letter_id()) {

        current = rule.pop_last_from_letter(current);
        assert(!current->is_nonterminal());

        if (next->is_valid() && current != next) {
          assert(current->last_terminal_letter_id() ==
                 next->first_terminal_letter_id());
          assert(next->is_nonterminal());
          next = rule.pop_first_from_letter(next);
        }

        assert(current->is_valid());

        assert(!next->is_valid() || current == next);
      } else {
        current = next;
      }
    }
  }
}


std::vector<LetterPosition> JezRules::list_blocks() {
  std::vector<LetterPosition> blocks;

  for (auto& rule : rules_) {
    for (auto letter = rule.begin(); letter != rule.end(); ++letter) {
      if (letter->is_power()) {
        blocks.push_back(LetterPosition(&rule, letter));
      }
    }
  }

  std::stable_sort(
    blocks.begin(),
    blocks.end(),
    [](const LetterPosition& first, const LetterPosition& second) -> bool {
      return first.letter_->terminal_id() < second.letter_->terminal_id() ||
          (first.letter_->terminal_id() == second.letter_->terminal_id() &&
           mad_sorts::reverse_bit_mpz_less(first.letter_->terminal_power(), second.letter_->terminal_power()));
    }
  );

  return blocks;
}

void JezRules::compress_blocks(const std::vector<LetterPosition>& blocks) {
  //LetterPower last_power = 0;
  TerminalId last_id = 0;

  Vertex current_terminal_vertex;
  std::vector<std::pair<LetterPower, Vertex*>> current_terminal_powers_;

  auto block = blocks.begin();
  while (block != blocks.end()) {
    assert(!block->letter_->is_nonterminal());
    if (!block->letter_->is_valid()) {
      ++block;
      continue;
    }

    assert(!block->rule_->empty());
    assert(block->letter_->is_power());

    if (last_id != block->letter_->terminal_id()) {
      last_id = block->letter_->terminal_id();
      current_terminal_vertex = terminal_vertices_[last_id];
      assert(current_terminal_powers_.empty());
    }

    if (current_terminal_powers_.empty() ||
          current_terminal_powers_.back().first !=
          block->letter_->terminal_power()
    ) {
      this->terminal_vertices_.insert(
        this->terminal_vertices_.end(),
        std::make_pair(
          next_fresh_terminal(),
          Vertex()
        )
      );

      current_terminal_powers_.push_back(
        std::make_pair(
          (block->letter_->terminal_power()),
          &(std::prev(this->terminal_vertices_.end())->second)
        )
      );
    }

    block->rule_->compress_power(
        block->letter_,
        last_terminal()
    );

    ++block;

    if (block == blocks.end() ||
        last_id != block->letter_->terminal_id()
    ) {
      bool continue_iterations = true;
      while (continue_iterations) {
        Vertex last_vertex;
        Vertex last_power = current_terminal_vertex;
        continue_iterations = false;
        for (auto& terminal_power : current_terminal_powers_) {
          if ((terminal_power.first & 1) != 0) {
            if (*terminal_power.second != last_vertex) {
              last_vertex = *terminal_power.second;
              last_power = NonterminalVertex(last_vertex, current_terminal_vertex);
            }
            *terminal_power.second = last_power;
          }

          terminal_power.first >>= 1;

          if (terminal_power.first != 0) {
            continue_iterations = true;
          }
        }
        current_terminal_vertex = NonterminalVertex(
          current_terminal_vertex,
          current_terminal_vertex
        );
      }
      current_terminal_powers_.clear();
    }
  }
}

void JezRules::empty_cleanup() {
  for (auto& rule : rules_) {
    for (auto current = rule.begin(); current != rule.end(); ) {
      auto next = std::next(current);
      if (current->is_empty_nonterminal()) {
        rule.remove_empty_letter(current);
      }
      current = next;
    }
  }
}

OneStepPairs::OneStepPairs(JezRules* rules)
  : rules_(rules)
{
  std::vector<std::tuple<
      TerminalId, //first letter
      TerminalId, //second letter
      LetterPosition>> all_pairs; //reference to position

  for (auto& rule : rules->rules_) {
    for (
        auto current = rule.begin(), next = std::next(rule.begin());
        next != rule.end();
        current = next, ++next
    ) {
      if (
          current->last_terminal_letter_id() != next->first_terminal_letter_id()
          //&& !current->is_power() && !next->is_power()
      ) {
        all_pairs.emplace_back(
            current->last_terminal_letter_id(),
            next->first_terminal_letter_id(),
            LetterPosition(&rule, current)
        );
      }
    }
  }

  std::stable_sort(
      all_pairs.begin(),
      all_pairs.end(),
      [] (const std::tuple<TerminalId, TerminalId, LetterPosition>& first,
          const std::tuple<TerminalId, TerminalId, LetterPosition>& second) {
        return std::get<0>(first) < std::get<0>(second) ||
            (std::get<0>(first) == std::get<0>(second) &&
             std::get<1>(first) < std::get<1>(second));
      }
  );

  if (all_pairs.empty()) {
    return;
  }

  auto left_letter_iterator = pairs_.end();
  auto right_letter_iterator = pairs_.end();

  decltype(right_letter_iterator->left_letters_.end()) left_list_current_letter;
  decltype(left_letter_iterator->right_letters_.end()) right_list_current_letter;

  for (const auto& pair: all_pairs) {
    while (
        left_letter_iterator != pairs_.end() &&
        left_letter_iterator->id_ < std::get<0>(pair)
    ) {
      ++left_letter_iterator;
      right_letter_iterator = pairs_.begin(); //TODO:optimize this
      right_list_current_letter = left_letter_iterator->right_letters_.begin();
      left_list_current_letter = right_letter_iterator->left_letters_.begin();
    }
    if (left_letter_iterator == pairs_.end() ||
        left_letter_iterator->id_ != std::get<0>(pair)
    ) {
      left_letter_iterator = pairs_.emplace(left_letter_iterator, std::get<0>(pair));
      right_letter_iterator = pairs_.begin();
      right_list_current_letter = left_letter_iterator->right_letters_.begin();
      left_list_current_letter = right_letter_iterator->left_letters_.begin();
    }

    while (
        right_list_current_letter != left_letter_iterator->right_letters_.end() &&
        right_list_current_letter->id_ < std::get<1>(pair)
    ) {
      ++right_list_current_letter;
    }

    if (
        right_list_current_letter == left_letter_iterator->right_letters_.end() ||
        right_list_current_letter->id_ != std::get<1>(pair)
    ) {
      right_list_current_letter = left_letter_iterator->right_letters_.insert(
          right_list_current_letter,
          std::get<1>(pair)
      );
    }

    right_list_current_letter->occurencies.push_back(std::get<2>(pair));

    while (
        right_letter_iterator != pairs_.end() &&
        right_letter_iterator->id_ < std::get<1>(pair)
    ) {
      ++right_letter_iterator;
      left_list_current_letter = right_letter_iterator->left_letters_.begin();
    }

    if (right_letter_iterator == pairs_.end() ||
        right_letter_iterator->id_ != std::get<1>(pair)
    ) {
      right_letter_iterator = pairs_.insert(
          right_letter_iterator,
          std::get<1>(pair)
      );

      left_list_current_letter = right_letter_iterator->left_letters_.begin();
    }

    while (
        left_list_current_letter != right_letter_iterator->left_letters_.end() &&
        *left_list_current_letter < std::get<0>(pair)
    ) {
      ++left_list_current_letter;
    }

    if (
        left_list_current_letter == right_letter_iterator->left_letters_.end() ||
        *left_list_current_letter != std::get<0>(pair)
    ) {
      left_list_current_letter = right_letter_iterator->left_letters_.insert(
          left_list_current_letter,
          std::get<0>(pair)
      );
    }
  }
}

std::tuple<std::vector<unsigned char>, std::vector<unsigned char>>
OneStepPairs::greedy_pairs() const {
  std::vector<unsigned char> left_letters;
  std::vector<unsigned char> right_letters;

  for (auto& letter_pairs : pairs_) {
    if (letter_pairs.right_letters_.empty() && letter_pairs.left_letters_.empty()) {
      continue;
    }

    size_t new_pairs_if_left = 0;

    for (const auto& letter : letter_pairs.right_letters_) {
      assert(letter.id_ != letter_pairs.id_);
      if (letter.id_ > letter_pairs.id_ || right_letters[letter.id_]) {
        ++new_pairs_if_left;
      }
    }

    size_t new_pairs_if_right = 0;
    for (const auto& letter : letter_pairs.left_letters_) {
      assert(letter != letter_pairs.id_);
      if (letter > letter_pairs.id_ || left_letters[letter]) {
        ++new_pairs_if_right;
      }
    }

    if (left_letters.empty()) {
      left_letters.resize(rules_->last_terminal() + 1, 0);
      right_letters.resize(rules_->last_terminal() + 1, 0);
    }

    if (new_pairs_if_left >= new_pairs_if_right) {
      left_letters[letter_pairs.id_] = 1;
    } else {
      right_letters[letter_pairs.id_] = 1;
    }
  }

  return std::make_tuple(std::move(left_letters), std::move(right_letters));
}

TerminalId OneStepPairs::compress_pair(
    TerminalId first,
    TerminalId second,
    const std::vector<LetterPosition>& occurencies
) {
  TerminalId pair_terminal_id = rules_->next_fresh_terminal();

  for (auto& occurence : occurencies) {
    auto first_letter = occurence.letter_;
    auto second_letter = std::next(first_letter);

    assert(first_letter->is_valid());
    assert(!occurence.rule_->empty());
    assert(!first_letter->is_nonterminal());
    assert(!first_letter->is_power());
    assert(second_letter != occurence.rule_->end());
    assert(!second_letter->is_nonterminal());
    assert(!second_letter->is_power());
    assert(first_letter->terminal_id() == first);
    assert(second_letter->terminal_id() == second);

    occurence.rule_->compress_pair(
        first_letter,
        second_letter,
        pair_terminal_id
    );
  }
  return pair_terminal_id;
}

void OneStepPairs::remove_crossing(
    const std::vector<unsigned char>& left_letters,
    const std::vector<unsigned char>& right_letters
) {
  std::vector<std::tuple<TerminalId, TerminalId, LetterPosition>> occurencies;

  for (auto & rule : rules_->rules_) {
    if (rule.empty()) {
      continue;
    }
    auto current = rule.begin();

    assert(!current->is_empty_nonterminal());

    auto next = std::next(current);
    while (current != rule.end() && next != rule.end()) {
      assert(!next->is_empty_nonterminal());

      assert(current->last_terminal_letter_id() >= 0);
      assert(next->first_terminal_letter_id() >= 0);


      if (current->last_terminal_letter_id() < left_letters.size() &&
          left_letters.at(current->last_terminal_letter_id()) &&
          next->first_terminal_letter_id() < right_letters.size() &&
          right_letters.at(next->first_terminal_letter_id()))
      {
        current = rule.pop_last_from_letter(current);

        assert(next->is_valid());
        next = rule.pop_first_from_letter(next);

        assert(current->is_valid());
        assert(!current->is_empty_nonterminal());
        assert(next->is_valid());
        assert(!next->is_empty_nonterminal());
        assert(current->last_terminal_letter_id() < left_letters.size() &&
               left_letters.at(current->last_terminal_letter_id()));
        assert(next->first_terminal_letter_id() < right_letters.size() &&
               right_letters.at(next->first_terminal_letter_id()));
        assert(std::next(current) == next);

        occurencies.emplace_back(current->last_terminal_letter_id(),
                                 next->first_terminal_letter_id(),
                                 LetterPosition(&rule, current));

      }
      current = next;
      next = std::next(current);
    }
  }

  std::stable_sort(
      occurencies.begin(),
      occurencies.end(),
      [] (const std::tuple<TerminalId, TerminalId, LetterPosition>& first,
          const std::tuple<TerminalId, TerminalId, LetterPosition>& second) {
        return std::get<0>(first) < std::get<0>(second) ||
            (std::get<0>(first) == std::get<0>(second) &&
             std::get<1>(first) < std::get<1>(second));
      }
  );

  auto current_pair = occurencies.begin();

  for (auto& pair : pairs_) {
    if (current_pair != occurencies.end() &&
        pair.id_ >= std::get<0>(*current_pair)) {
      assert(pair.id_ == std::get<0>(*current_pair));
      for (auto & right_letter : pair.right_letters_) {
        right_letter.occurencies.clear();
        if (current_pair != occurencies.end() &&
            pair.id_ >= std::get<0>(*current_pair) &&
            right_letter.id_ >= std::get<1>(*current_pair)) {
          assert(right_letter.id_ == std::get<1>(*current_pair));
          while (pair.id_ == std::get<0>(*current_pair) &&
              right_letter.id_ == std::get<1>(*current_pair)) {
            right_letter.occurencies.push_back(std::get<2>(*current_pair));
            ++current_pair;

            assert(!right_letter.occurencies.back().rule_->empty());
            assert(!right_letter.occurencies.back().letter_->is_nonterminal());
            assert(!right_letter.occurencies.back().letter_->is_power());
            assert(std::next(right_letter.occurencies.back().letter_) !=
                right_letter.occurencies.back().rule_->begin());
            assert(!std::next(right_letter.occurencies.back().letter_)->is_nonterminal());
            assert(!std::next(right_letter.occurencies.back().letter_)->is_power());

            if (current_pair == occurencies.end()) {
              break;
            }
          }
        }
      }
    } else {
      for (auto & right_letter : pair.right_letters_) {
        right_letter.occurencies.clear();
      }
    }
  }
}


void OneStepPairs::compress_pairs_from_letter_lists(
    const std::vector<unsigned char>& left_letters,
    const std::vector<unsigned char>& right_letters
) {
  for (auto& pair : pairs_) {
    if (pair.id_ < left_letters.size() && left_letters.at(pair.id_)) {
      auto right_letter = pair.right_letters_.begin();

      while (right_letter != pair.right_letters_.end()) {
        if (right_letter->id_ < right_letters.size() && right_letters.at(right_letter->id_)) {
          auto terminal = compress_pair(
              pair.id_,
              right_letter->id_,
              right_letter->occurencies
          );

          rules_->terminal_vertices_.insert(std::make_pair(
              terminal,
              NonterminalVertex(
                  rules_->terminal_vertices_[pair.id_],
                  rules_->terminal_vertices_[right_letter->id_]
              )
          ));

          right_letter = pair.right_letters_.erase(right_letter);
        } else {
          ++right_letter;
        }
      }
    }

    if (pair.id_ < right_letters.size() && right_letters.at(pair.id_)) {
      auto left_letter = pair.left_letters_.begin();

      while (left_letter != pair.left_letters_.end()) {
        if (*left_letter < left_letters.size() && left_letters.at(*left_letter)) {
          left_letter = pair.left_letters_.erase(left_letter);
        } else {
          ++left_letter;
        }
      }
    }
  }
}

Vertex JezRules::get_exponent(Vertex vertex, LetterPower power) {
  Vertex result;

  assert(power > 0);

  while (power != 0) {
    if ((power & 1) != 0) {
      if (result) {
        result = NonterminalVertex(result, vertex);
      } else {
        result = vertex;
      }
    }
    vertex = NonterminalVertex(vertex, vertex);
    power >>= 1;
  }

  return result;
}

void RuleLetter::debug_print(std::ostream* os) const {
  if (this->is_nonterminal()) {
    (*os) << "(" << this->nonterminal_rule()->debug_id << ")";
  } else {
    (*os) << this->terminal_id() << "^" << this->terminal_power();
  }
}

void Rule::debug_print(::std::ostream* os) const {
  for (auto& letter : letters_) {
    letter.debug_print(os);
  }
}

void Rule::debug_print_exposed(::std::ostream* os) const {
  for (auto& letter : letters_) {
    if (letter.is_nonterminal()) {
      letter.nonterminal_rule_->debug_print_exposed(os);
    } else {
      for (LetterPower i = 0; i < letter.terminal_power(); ++i) {
        (*os) << letter.terminal_id() << ',';
      }
    }
  }
}

void JezRules::debug_print(std::ostream* os) const {
  for (auto& rule : rules_) {
    (*os) << rule.debug_id << ": ";
    (*os) << rule.first_terminal_id() << ".." <<
        rule.last_terminal_id() << ": ";
    for (auto& letter : rule) {
      letter.debug_print(os);
      (*os) << ", ";
    }
    (*os) << std::endl;
  }
}

Vertex normal_form(Vertex root) {
  if (root.height() < 2) {
    return root;
  }
  JezRules rules(root);
  Rule& root_rule = *(rules.vertex_rules_[root]);

  while (root_rule.size() > 1 || (
           root_rule.size() == 1 && (
             root_rule.begin()->is_power() ||
             root_rule.begin()->is_nonterminal()
           )
        )) {

//    std::cout << "\n=================\n\nCurrent rules:" << std::endl;
//    rules.debug_print(&std::cout);
    OneStepPairs pairs(&rules);

    rules.remove_crossing_blocks();

//
//    std::cout << "Rules after RemCrBlocks: " << std::endl;
//    rules.debug_print(&std::cout);

    auto blocks = rules.list_blocks();
//    std::cout << "\nFound blocks: " << std::endl;
//    for (auto& block : blocks) {
//      std::cout << block.rule_->debug_id << ':';
//      block.letter_->debug_print(&std::cout);
//      std::cout << std::endl;
//    }

    rules.compress_blocks(blocks);

//    std::cout << "Rules after CompressBlocks: " << std::endl;
//    rules.debug_print(&std::cout);

    std::vector<unsigned char> left_letters, right_letters;

    std::tie(left_letters, right_letters) = pairs.greedy_pairs();
//    std::cout << "\nGreedyPairs:\nLeft: ";
//
//    for (auto& terminal : left_letters) {
//      std::cout << terminal << ',';
//    }
//    std::cout << "\nRight: ";
//    for (auto& terminal : right_letters) {
//      std::cout << terminal << ',';
//    }
//    std::cout << std::endl;

    while (!left_letters.empty()) {
      pairs.remove_crossing(left_letters, right_letters);
      pairs.compress_pairs_from_letter_lists(left_letters, right_letters);
//      std::cout << "Rules after first compression: " << std::endl;
//      rules.debug_print(&std::cout);
      pairs.remove_crossing(right_letters, left_letters);
      pairs.compress_pairs_from_letter_lists(right_letters, left_letters);
//      std::cout << "Rules after second compression: " << std::endl;
//      rules.debug_print(&std::cout);
      std::tie(left_letters, right_letters) = pairs.greedy_pairs();
//      std::cout << "\nGreedyPairs:\nLeft: ";
//
//      for (auto& terminal : left_letters) {
//        std::cout << terminal << ',';
//      }
//      std::cout << "\nRight: ";
//      for (auto& terminal : right_letters) {
//        std::cout << terminal << ',';
//      }
//      std::cout << std::endl;
    }

    rules.empty_cleanup();
  }

  Rule::collect_garbage();
  return rules.terminal_vertices_[root_rule.first_terminal_id()];
}

}
}
}



