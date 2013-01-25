#include "SLPSet.h"
#include <gmp.h>

namespace crag {

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
  PathState current = current_path_.back();
  current_path_.pop_back();

  if (inspection_ended()) {
    return;
  }

  const PathState& parent = current_path_.back();

  if ( current.is_right_child || !parent.vertex.has_right_child()) {
    //We have already visited all right children of the parent vertex, stop on parent
  } else {
    //We have visited all left children of parent, going to the right
    current_path_.push_back(PathState({parent.vertex.right_child(), true}));
    goto_leftmost_terminal();
  }
}

void SLPPostorderInspector::goto_leftmost_terminal() {
  while (current_path_.back().vertex != SLPVertex::Null) {
    const PathState & current = current_path_.back();
    if (current.vertex.has_left_child()) {
      current_path_.push_back(PathState({current.vertex.left_child(), false}));
    } else {
      //Don't have any left children, then go right
      current_path_.push_back(PathState({current.vertex.right_child(), true}));
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

  if (pattern.length() == 0 || text.length() == 0) {
    return NO_MATCHES;
  }

  auto match_result_iterator = match_table_.find(std::make_pair(pattern, text));

  if (match_result_iterator != match_table_.end()) { //if already calculated
    return match_result_iterator->second;
  }

  if (pattern.length() > text.length()) {
    return NO_MATCHES;
  }

  MatchResultSequence match_result;

  if (pattern == text) { //If we are checking the same vertex
    match_result = {0, 1, 1};
  } else if (pattern.length() == 1) {//Trivial case
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
    if (pattern.is_negative()) {//minor optimization
      auto inversed_result = this->matches(pattern.negate(), text.negate());

      if (inversed_result.count <= 0) {
        match_result = NO_MATCHES;
      } else {
        match_result = {
            text.length() - pattern.length() - inversed_result.start - inversed_result.step * (inversed_result.count - 1),
            inversed_result.step,
            inversed_result.count
        };
      }
    } else if (pattern.left_child().length() >= pattern.right_child().length()) {//Right child is smaller
      match_result = internal::nontrivial_match(
          pattern.left_child(),
          pattern.right_child(),
          true,
          text,
          this
      );
    } else {
      match_result = internal::nontrivial_match(
          pattern.right_child(),
          pattern.left_child(),
          false,
          text,
          this
      );
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

  SLPMatchingInspector large_part_hunter(
      large_pattern_part.length(), //pattern length
      text, //text
      large_pattern_part_left_bound, //the leftmost letter of text to consider
      large_pattern_part.length() * 2 + small_pattern_part.length()
  );

  SLPMatchingTable::MatchResultSequence result;

  while(!large_part_hunter.inspection_ended()) {
    auto large_part_matches = internal::local_search(large_pattern_part, &large_part_hunter, matching_table);

    if (large_part_matches.count > 0) {
      SLPMatchingTable::MatchResultSequence small_part_candidates;
      if (small_pattern_is_after) {
        small_part_candidates = {
            large_part_matches.start + large_pattern_part.length(),
            large_part_matches.step,
            large_part_matches.count
        };
      } else {
        small_part_candidates = {
            large_part_matches.start - small_pattern_part.length(),
            large_part_matches.step,
            large_part_matches.count
        };

        if (small_part_candidates.start < 0) {
          LongInteger count_left_cut;
          LongInteger new_start;
          mpz_fdiv_qr(count_left_cut.get_mpz_t(), new_start.get_mpz_t(), small_part_candidates.start.get_mpz_t(), small_part_candidates.step.get_mpz_t());
          small_part_candidates.start = new_start;
          small_part_candidates.count += count_left_cut;
        }
      }


      const LongInteger& first_candidate_start = small_part_candidates.start;
      LongInteger last_candidate_start = small_part_candidates.start + small_part_candidates.step * (small_part_candidates.count - 1);

      LongInteger seaside_candidates_bound = (small_pattern_is_after ? last_candidate_start - small_pattern_part.length(): first_candidate_start);

      SLPMatchingInspector seaside_hunter(
          small_pattern_part.length(),
          text,
          seaside_candidates_bound,
          2 * small_pattern_part.length()
      );

      auto seaside_matches = internal::local_search(small_pattern_part, &seaside_hunter, matching_table);
      auto seaside_approved_candidates = internal::intersect_sequences(small_part_candidates, seaside_matches);

      LongInteger continental_candidate_start = small_pattern_is_after ? first_candidate_start : last_candidate_start;
      SLPMatchingInspector continental_hunter(
          small_pattern_part.length(),
          text,
          continental_candidate_start,
          small_pattern_part.length()
      );

      auto continental_matches = internal::local_search(small_pattern_part, &continental_hunter, matching_table);
      SLPMatchingTable::MatchResultSequence approved_candidates;
      if (continental_matches.count > 0) {
        SLPMatchingTable::MatchResultSequence continental_approved_candidates = small_part_candidates;
        LongInteger small_pattern_part_length_in_large_steps;
        LongInteger small_pattern_shift;
        mpz_cdiv_qr(small_pattern_part_length_in_large_steps.get_mpz_t(), small_pattern_shift.get_mpz_t(), small_pattern_part.length().get_mpz_t(), large_part_matches.step.get_mpz_t());
        continental_approved_candidates.count -= small_pattern_part_length_in_large_steps;
        if (!small_pattern_is_after) {
          continental_approved_candidates.start += small_pattern_part.length() - small_pattern_shift;
        }
        approved_candidates = internal::join_sequences(continental_approved_candidates, seaside_approved_candidates);
      } else {
        approved_candidates = seaside_approved_candidates;
      }

      if (approved_candidates.count > 0 && small_pattern_is_after) {
        approved_candidates.start -= large_pattern_part.length();
        if (approved_candidates.start < 0) {
          LongInteger count_left_cut;
          LongInteger new_start;
          mpz_fdiv_qr(count_left_cut.get_mpz_t(), new_start.get_mpz_t(), approved_candidates.start.get_mpz_t(), approved_candidates.step.get_mpz_t());
          approved_candidates.start = new_start;
          approved_candidates.count += count_left_cut;
        }
      }

      result = internal::join_sequences(result, approved_candidates);
    }
  }

  if (result.count == 0) {
    return SLPMatchingTable::NO_MATCHES;
  }

  if (result.count == 1) {
    result.step = 1;
  }

  return result;
}


void internal::SLPMatchingInspector::go_further() {
    while (position_is_suitable()) {
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
  if (first.step == 0 || first.count <= 0) {
    return second;
  }
  if (second.step == 0  || second.count <= 0) {
    return first;
  }
  //if one of the sequences is just one element, then the behavior should be different
  if (first.count == 1) {
    std::swap(first, second);
  }

  if (second.count == 1) {
    if (first.count == 1) {
      if (first.start < second.start) {
        first.step = second.start - first.start;
        first.count = 2;
        return first;
      } else if (second.start < first.start) {
        second.step = first.start - second.start;
        second.count = 2;
        return second;
      } else {
        second.step = 1;
        second.count = 1;
        return second;
      }
    } else {
      LongInteger distance_between_starts = second.start - first.start;
      LongInteger steps_between_starts;
      LongInteger distance_badness;

      mpz_fdiv_qr(steps_between_starts.get_mpz_t(), distance_badness.get_mpz_t(), distance_between_starts.get_mpz_t(), first.step.get_mpz_t());

      if (distance_badness != 0) {
        return SLPMatchingTable::NO_MATCHES;
      }

      if (steps_between_starts == -1) {
        first.start -= first.step;
        first.count += 1;

        return first;
      } else if (steps_between_starts == first.count) {
        first.count += 1;

        return first;
      } else if (steps_between_starts >= 0 && steps_between_starts < first.count) {
        return first;
      } else {
        return SLPMatchingTable::NO_MATCHES;
      }
    }
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
    second.count += steps_inside_starts_interval;
    second.start = first.start;
    return second;
  } else {
    return first;
  }
}

SLPMatchingTable::MatchResultSequence internal::intersect_sequences(SLPMatchingTable::MatchResultSequence first,
                                                     SLPMatchingTable::MatchResultSequence second) {
  if (first.step == 0 || first.count == 0 || second.step == 0 || second.count == 0) {
    return SLPMatchingTable::NO_MATCHES;
  }

  if (first.start > second.start) {
    std::swap(first, second);
  }

  LongInteger first_step_coefficient;
  LongInteger second_step_coefficient;
  LongInteger steps_gcd;

  mpz_gcdext(steps_gcd.get_mpz_t(), first_step_coefficient.get_mpz_t(), second_step_coefficient.get_mpz_t(), first.step.get_mpz_t(), second.step.get_mpz_t());

  LongInteger starts_difference = first.start - second.start;

  if (!mpz_divisible_p(starts_difference.get_mpz_t(), steps_gcd.get_mpz_t())) {
    return SLPMatchingTable::NO_MATCHES;
  }

  SLPMatchingTable::MatchResultSequence result;

  result.step = (first.step * second.step / steps_gcd);

  LongInteger remainder;
  second_step_coefficient *= starts_difference;
  mpz_fdiv_r(remainder.get_mpz_t(), second_step_coefficient.get_mpz_t(), first.step.get_mpz_t());

  result.start = second.start  + remainder * second.step / steps_gcd;

  LongInteger last_first = first.start + first.step * (first.count - 1);
  LongInteger last_second = second.start + second.step * (second.count - 1);

  if (last_first > last_second) {
    last_first.swap(last_second);
  }

  last_first -= result.start;

  mpz_fdiv_q(result.count.get_mpz_t(), last_first.get_mpz_t(), result.step.get_mpz_t());
  result.count += 1;

  if (result.count <= 0) {
    return SLPMatchingTable::NO_MATCHES;
  }

  if (result.count == 1) {
    result.step = 1;
  }

  return result;
}

SLPMatchingTable::MatchResultSequence internal::local_search(const SLPVertex& pattern, SLPMatchingInspector* inspector, SLPMatchingTable* matching_table) {
  SLPMatchingTable::MatchResultSequence current_result = SLPMatchingTable::NO_MATCHES;
  while (!inspector->inspection_ended()) {
    auto current_match = matching_table->matches(pattern, inspector->current_text());
    //Adjust match start
    current_match.start += (inspector->current_distance_from_left_border() + inspector->left_border());
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
    LongInteger current_match_last_element_end = current_match.start + (current_match.count - 1) * current_match.step + pattern.length();
    if (current_match.count > 0 && current_match_last_element_end > inspector->left_border() + inspector->interval_length()) {
      LongInteger count_reduce;
      LongInteger current_end_offset = current_match_last_element_end - inspector->left_border() - inspector->interval_length();
      mpz_cdiv_q(count_reduce.get_mpz_t(), current_end_offset.get_mpz_t(), current_match.step.get_mpz_t());
      current_match.count -= count_reduce;
    }

    auto joined_result = join_sequences(current_result, current_match);
    if (current_result.count > 0 &&
        (joined_result.count <= 0 ||
         joined_result.step > pattern.length() //any two consequent matches must have some common point
        )) {
      return current_result;
    }
    current_result = std::move(joined_result);
    inspector->go_next();
  }
  return current_result;
}

::std::ostream& operator<<(::std::ostream& os, const SLPMatchingTable::MatchResultSequence& match) {
  return os << '{' << match.start << ',' << match.step << ',' << match.count << '}';  // whatever needed to print bar to os
}
} //namespace crag

void std::swap(crag::SLPMatchingTable::MatchResultSequence& first, crag::SLPMatchingTable::MatchResultSequence& second) {
  mpz_swap(first.count.get_mpz_t(), second.count.get_mpz_t());
  mpz_swap(first.start.get_mpz_t(), second.start.get_mpz_t());
  mpz_swap(first.step.get_mpz_t(), second.step.get_mpz_t());
}
