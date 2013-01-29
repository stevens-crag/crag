#include "SLPSet.h"

namespace crag {

const SLPVertex SLPVertex::Null; //define static in cpp file

inline BasicVertex::BasicVertex()
  : left_child(SLPVertex::Null),
    right_child(SLPVertex::Null),
    terminal_symbol(INVALID_TERMINAL),
    length(0),
    height(0)
{ }

SLPVertex SLPVertex::terminal_vertex(const TerminalSymbol& terminal_symbol) {
  SLPVertex vertex(std::make_shared<BasicVertex>(), false);

  vertex.ptr_->terminal_symbol = terminal_symbol;
  vertex.ptr_->length = 1;
  vertex.ptr_->height = 1;

  return vertex;
}

SLPVertex SLPVertex::concatenate(const SLPVertex& left_child, const SLPVertex& right_child) {
  SLPVertex vertex(std::make_shared<BasicVertex>(), false);

  vertex.ptr_->left_child = left_child;
  vertex.ptr_->right_child = right_child;
  vertex.ptr_->length = left_child.length() + right_child.length();
  vertex.ptr_->height = std::max(left_child.height(), right_child.height()) + 1;

  return vertex;
}

void SLPPostorderInspector::go_to_next_vertex() {
  //TODO: maybe raise an exception if current_path is empty already?
  const SLPVertex current = current_path_.back();
  current_path_.pop_back();

  if (inspection_ended()) {
    return;
  }

  const SLPVertex& parent = current_path_.back();

  if ( parent.right_child() == current || !parent.has_right_child()) {
    //We have already visited all right children of the parent vertex, stop on parent
  } else {
    //We have visited all left children of parent, going to the right
    current_path_.push_back(parent.right_child());
    goto_leftmost_terminal();
  }
}

void SLPPostorderInspector::goto_leftmost_terminal() {
  while (current_path_.back() != SLPVertex::Null) {
    const SLPVertex & current = current_path_.back();
    if (current.has_left_child()) {
      current_path_.push_back(current.left_child());
    } else {
      //Don't have any left children, then go right
      current_path_.push_back(current.right_child());
    }
  }

  current_path_.pop_back();
}

SLPProducedWordIterator& SLPProducedWordIterator::operator ++() {
  ++length_;
  inspector_.go_to_next_vertex();

  while (!inspector_.inspection_ended() && !inspector_.current_vertex().is_terminal()) {
    inspector_.go_to_next_vertex();
  }

  return *this;
}

SLPProducedWord::value_type SLPProducedWord::operator [](LongInteger index) const {
  SLPVertex current_vertex = root_;

  while (!current_vertex.is_terminal()) {
    if (index < current_vertex.left_child().length()) { //Go left
      current_vertex = current_vertex.left_child();
    } else {
      index -= current_vertex.left_child().length();
      current_vertex = current_vertex.right_child();
    }
  }

  return current_vertex;
}

const SLPMatchingTable::MatchResultSequence SLPMatchingTable::NO_MATCHES;

bool terminals_equal(const SLPVertex& lhs, const SLPVertex& rhs) {
  return lhs.terminal_symbol() == rhs.terminal_symbol() && lhs.is_negative() == rhs.is_negative();
}

SLPMatchingTable::MatchResultSequence SLPMatchingTable::matches(const SLPVertex& pattern,
                                              const SLPVertex& text) {
  auto match_result_iterator = match_table_.find(std::make_pair(pattern, text));

  if (match_result_iterator != match_table_.end()) { //if already calculated
    return match_result_iterator->second;
  }

  if (pattern.length() > text.length()) {
    return NO_MATCHES;
  }

  MatchResultSequence match_result;

  if (pattern.length() == 1) {//Trivial case
    SLPVertex pattern_vertex = SLPProducedWord(pattern)[0];
    if (text.length() == 1) {
      if (terminals_equal(pattern_vertex, SLPProducedWord(text)[0])) {
        match_result.count = 1;
        match_result.start = 0;
        match_result.step = 1;
      } else {
        match_result = NO_MATCHES;
      }
    } else {
      SLPProducedWord left_text_part = SLPProducedWord(text.left_child());
      SLPProducedWord right_text_part = SLPProducedWord(text.right_child());

      bool left_match = left_text_part.size() > 0 &&
          terminals_equal(pattern_vertex, left_text_part[left_text_part.size() - 1]);
      bool right_match = right_text_part.size() > 0 &&
          terminals_equal(pattern_vertex, right_text_part[0]);

      if (left_match) {
        match_result.count = right_match ? 2 : 1;
        match_result.start = left_text_part.size() - 1;
        match_result.step = 1;
      } else if (right_match) {
        match_result.count = 1;
        match_result.start = left_text_part.size();
        match_result.step = 1;
      } else {
        match_result = NO_MATCHES;
      }
    }
  } else  {//we have pattern.length > 1 => text.length > 1
    if (pattern.left_child().length() >= pattern.right_child().length()) {//Right child is smaller
      const SLPVertex& large_pattern_part = pattern.left_child();
      const SLPVertex& small_pattern_part = pattern.right_child();
      bool small_pattern_is_after = true;
    }
  }

  match_table_.insert(std::make_pair(std::make_pair(pattern, text), match_result));

  return match_result;
}

SLPMatchingTable::MatchResultSequence internal::nontrivial_match(const SLPVertex& large_pattern_part,
                                                       const SLPVertex& small_pattern_part,
                                                       bool small_pattern_is_after,
                                                       const SLPVertex& text,
                                                       SLPMatchingTable* matching_table) {
  LongInteger large_pattern_part_left_bound = text.split_point() - large_pattern_part.length();
  if (small_pattern_is_after) {
    large_pattern_part_left_bound -= small_pattern_part.length();
  }

  LongInteger large_pattern_part_right_bound = text.split_point() + large_pattern_part.length();

  if (!small_pattern_is_after) {
    large_pattern_part_right_bound += small_pattern_part.length();
  }

  SLPMatchingInspector large_part_hunter(
      large_pattern_part.length(), //pattern length
      text, //text
      large_pattern_part_left_bound, //the leftmost letter of text to consider
      large_pattern_part_right_bound //the rightmost letter of text to consider
  );

  SLPMatchingTable::MatchResultSequence result;

  while(!large_part_hunter.inspection_ended()) {
    auto large_part_matches = local_search(large_pattern_part, &large_part_hunter, matching_table);
    if (large_part_matches.count > 0) {
      if (small_pattern_is_after) {
        LongInteger large_pattern_first_match_end = large_part_matches.start + large_pattern_part.length();
        LongInteger large_pattern_last_match_end = large_pattern_first_match_end + large_part_matches.count * large_part_matches.step;

        SLPMatchingInspector seaside_hunter(
            small_pattern_part.length(),
            text,
            large_pattern_last_match_end - small_pattern_part.length(),
            large_pattern_last_match_end + small_pattern_part.length()
        );

        auto seaside_matches = local_search(small_pattern_part, &seaside_hunter, matching_table);

        SLPMatchingInspector continental_hunter(
            small_pattern_part.length(),
            text,
            large_pattern_first_match_end,
            large_pattern_first_match_end + small_pattern_part.length()
        );

        auto continental_matches = local_search(small_pattern_part, &continental_hunter, matching_table);

        //TODO: join all matches


      }
    }
  }

  return result;
}


void internal::SLPMatchingInspector::go_further() {
    while (position_is_suitable(text_position_)) {
      parent_stack_.push_back(text_position_);
      text_position_.vertex = text_position_.vertex.left_child();
    }
    if (!parent_stack_.empty()) {
      text_position_ = parent_stack_.back();
      parent_stack_.pop_back();
    } else {
      text_position_ = {SLPVertex::Null, 0};
    }
}

SLPMatchingTable::MatchResultSequence internal::join_sequences(SLPMatchingTable::MatchResultSequence first, SLPMatchingTable::MatchResultSequence second) {
  if (first.step == 0) {
    return second;
  }
  if (second.step == 0) {
    return first;
  }
  if (first.step != second.step) {
    return SLPMatchingTable::NO_MATCHES;
  }

  if (first.start > second.start) { //we assume that the first sequence begin before the second
    std::swap(first, second);
  }

  LongInteger distance_between_starts = second.start - first.start;
  LongInteger steps_inside_starts_interval;
  LongInteger residue; //Should be 0

  mpz_fdiv_qr(steps_inside_starts_interval.get_mpz_t(), residue.get_mpz_t(), distance_between_starts.get_mpz_t(), first.step.get_mpz_t());

  if (residue != 0) { //starts are not coherent with step
    return SLPMatchingTable::NO_MATCHES;
  }

  if (first.count < steps_inside_starts_interval) { //first sequence ends before second starts
    return SLPMatchingTable::NO_MATCHES;
  }

  if (first.count < second.count + steps_inside_starts_interval) {//Second sequence ends after first
    first.count = second.count + steps_inside_starts_interval;
    return first;
  }

  return first;
}


SLPMatchingTable::MatchResultSequence internal::local_search(const SLPVertex& pattern, SLPMatchingInspector* inspector, SLPMatchingTable* matching_table) {
  SLPMatchingTable::MatchResultSequence current_result = SLPMatchingTable::NO_MATCHES;
  while (!inspector->inspection_ended()) {
    auto current_match = matching_table->matches(pattern, inspector->current_text());
    //cut the sequence to begin after the large_pattern_part_left_bound
    if (current_match.count > 0 && current_match.start < inspector->left_border()) {
      LongInteger count_reduce; //The number of steps which are outside of the boundary
      LongInteger correct_begin_offset; //Distance between boundary and first match
      LongInteger current_begin_offset = inspector->left_border() - current_match.start;
      mpz_cdiv_qr(count_reduce.get_mpz_t(), correct_begin_offset.get_mpz_t(), current_begin_offset.get_mpz_t(), current_match.step.get_mpz_t());
      current_match.count -= count_reduce;
      current_match.start += current_begin_offset - correct_begin_offset; //They have different signs, because they are located on the different sides from left bound
    }
    //cut the sequence to end before large_pattern_part_right_bound
    LongInteger current_match_last_element_end = current_match.start + current_match.count * current_match.step + pattern.length();
    if (current_match.count > 0 && current_match_last_element_end > inspector->right_border()) {
      LongInteger count_reduce;
      LongInteger current_end_offset = current_match_last_element_end - inspector->right_border();
      mpz_cdiv_q(count_reduce.get_mpz_t(), current_end_offset.get_mpz_t(), current_match.step.get_mpz_t());
      current_match.count -= count_reduce;
    }
    auto joined_result = join_sequences(current_result, current_match);
    if (current_result != SLPMatchingTable::NO_MATCHES && joined_result == SLPMatchingTable::NO_MATCHES) {
      return current_result;
    }
    current_result = std::move(joined_result);
    inspector->go_next();
  }
  return current_result;
}


} //namespace crag
