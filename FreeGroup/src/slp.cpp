/*
 * slp.cpp
 *
 *  Created on: Feb 20, 2013
 *      Author: dpantele
 */

#include <iostream>

#include "slp.h"

namespace crag {
namespace slp {

const Vertex Vertex::Null;
const LongInteger Vertex::LongZero;

Vertex internal::BasicNonterminalVertex::negate() const {
  return NonterminalVertex(::std::make_shared<BasicNonterminalVertex>(
      ::std::shared_ptr<NonterminalVertexNodeData>(node_data_ptr_),
      !negate_node_
  ));
}

const ::std::hash<std::shared_ptr<internal::NonterminalVertexNodeData>> internal::BasicNonterminalVertex::ptr_hash;
size_t internal::BasicNonterminalVertex::last_vertex_id_;

const inspector::InspectorTask inspector::InspectorTask::DO_NOTHING = inspector::InspectorTask();

}

::std::ostream& operator<<(::std::ostream& out, const FiniteArithmeticSequence& sequence) {
  return out << sequence.first() << ":" << sequence.last() << ".." << sequence.step();
}

FiniteArithmeticSequence& FiniteArithmeticSequence::fit_into(const LongInteger& left_bound, const LongInteger& right_bound) {
  if (!*this) {
    return *this;
  }

  if (left_bound > last_ || right_bound < first_) {
    return *this = FiniteArithmeticSequence::Null;
  }

  //TODO: add static to speedup

  if (first_ < left_bound) {
    LongInteger current_begin_offset = left_bound - first_;
    mpz_cdiv_r(first_.get_mpz_t(), current_begin_offset.get_mpz_t(), step_.get_mpz_t());
    mpz_neg(first_.get_mpz_t(), first_.get_mpz_t());
    first_ += left_bound;
  }

  if (last_ > right_bound) {
    LongInteger current_end_offset = last_ - right_bound;
    mpz_cdiv_r(last_.get_mpz_t(), current_end_offset.get_mpz_t(), step_.get_mpz_t());
    last_ += right_bound;
  }

  if (first_ > last_) {
    *this = FiniteArithmeticSequence::Null;
  } else if (first_ == last_) {
    step_ = 1;
  }

  return *this;
}


FiniteArithmeticSequence& FiniteArithmeticSequence::join_with(const FiniteArithmeticSequence& other) {
  if (!*this) {
    return *this = other;
  }

  if (!other) {
    return *this;
  }

  if (this->step_ != other.step_ && other.last_ != other.first_ && this->last_ != this->first_) {
    return *this = FiniteArithmeticSequence();
  }

  if (first_ > other.first_) {
    return *this = ::std::move(FiniteArithmeticSequence(other).join_with(*this));
  }

  //Let make this a bit messy
  static mpz_t distance_between_starts,
         steps_inside_starts_interval,
         residue;  //Should be 0
  //mpz_inits(distance_between_starts, steps_inside_starts_interval, residue, NULL);

  mpz_sub(distance_between_starts, other.first_.get_mpz_t(), this->first_.get_mpz_t());
  if (this->first_ == this->last_) {
    if (other.first_ == other.last_) {
      if (mpz_cmp_ui(distance_between_starts, 0) != 0) {
        mpz_set(this->step_.get_mpz_t(), distance_between_starts);
      }
    } else {
      this->step_ = other.step_;
    }
  }

  mpz_fdiv_qr(steps_inside_starts_interval, residue, distance_between_starts, step_.get_mpz_t());

  if (mpz_cmp_ui(residue, 0) != 0) { //starts are not coherent with step
    *this = FiniteArithmeticSequence();
  } else if (other.first_ > this->last_ + this->step_) { //first sequence ends before second starts
    *this = FiniteArithmeticSequence();
  } else {
    if (other.last_ > this->last_) {
      this->last_ = other.last_;
    }
  }

  //mpz_clears(distance_between_starts, steps_inside_starts_interval, residue, NULL);
  return *this;
}

FiniteArithmeticSequence& FiniteArithmeticSequence::intersect_with(const FiniteArithmeticSequence& other) {
  if (!*this || !other) {
    return *this = FiniteArithmeticSequence();
  }

  if (other.first_ > this->first_) { //Ensure that this starts after other
    return *this = ::std::move(FiniteArithmeticSequence(other).intersect_with(*this));
  }
  //from this point we assign index 1 to the other sequence and index 2 to this sequence, so first_1 <= first_2

  static mpz_t temp, //static here gives 5% performance
        step_2_inverse,
        steps_gcd;
  //mpz_inits(temp, step_2_inverse, steps_gcd, intersection_end, NULL);

  //Here temp is just coefficient in front of other.step in extended gcd, not used
  mpz_gcdext(steps_gcd, temp, step_2_inverse, other.step_.get_mpz_t(), this->step_.get_mpz_t()); //gcd = u * step_1 + v * step_2
  //we know that the step of the result must be step_1 * step_2 / gcd.
  //the problem is to find minimal k >= 0 such that
  //first_2 + k * step_2 = first_1 (mod step_1)
  //then k = (first_1 - first_2) / gcd * (step_2/gcd)^{-1} (mod step_1/gcd)
  //note that step_2_inverse = (step_2/gcd)^{-1} (mod step_1/gcd)
  //so we just have to calculate (first_1 - first_2) * step_2_inverse % step_1 == k * gcd
  //after that the first element of result is first_2 + k * gcd * (step_2 / gcd)

  //But first, check if first_1 - first_2 is divisible by gcd
  //now temp is first_1 - first_2
  mpz_sub(temp, other.first_.get_mpz_t(), this->first_.get_mpz_t());

  if (!mpz_divisible_p(temp, steps_gcd)) { //if (first_2 - first_1) is not divisible by gcd, then sequences has completely different elements, not intersecting ever
    *this = FiniteArithmeticSequence();
  } else {
    mpz_mul(step_2_inverse, step_2_inverse, temp); //step_2_inverse = (first_1 - first_2) * step_2_inverse

    //we store step_2 / gcd in this->step_
    mpz_divexact(this->step_.get_mpz_t(), this->step_.get_mpz_t(), steps_gcd);

    //now we have to calculate the start of resulting sequence
    mpz_fdiv_r(temp, step_2_inverse, other.step_.get_mpz_t()); //temp = k * gcd = step_2_inverse % step_1
    mpz_addmul(this->first_.get_mpz_t(), temp, this->step_.get_mpz_t()); //first_result = first_2 + (k * gcd) * (step_2 / gcd)

    //now calculate the result step as (step_2 / gcd) * step_1
    this->step_ *= other.step_;

    if (this->last_ > other.last_) {
      this->last_ = other.last_;
    }

    mpz_sub(temp, this->last_.get_mpz_t(), this->first_.get_mpz_t());
    mpz_fdiv_r(temp, temp, this->step_.get_mpz_t());

    mpz_sub(this->last_.get_mpz_t(), this->last_.get_mpz_t(), temp);

    if (this->last_ < this->first_) {
      *this = FiniteArithmeticSequence();
    } else if (this->last_ == this->first_) {
      this->step_ = 1;
    }
  }

  //mpz_clears(temp, step_2_inverse, steps_gcd, intersection_end, NULL);

  return *this;
}

const FiniteArithmeticSequence FiniteArithmeticSequence::Null = FiniteArithmeticSequence();

namespace slp {

FiniteArithmeticSequence PatternMatchesGenerator::next_match() {
  FiniteArithmeticSequence result = FiniteArithmeticSequence::Null;

  while (!text_inspector_.stopped()) {
    //::std::cout << "Generator on ";
    //PrintTo(text_inspector_.vertex(), &::std::cout);
    //::std::cout << ::std::endl;
    FiniteArithmeticSequence new_match = matching_table_->matches(pattern_, text_inspector_.vertex());
    //::std::cout << "Match is " << new_match << ::std::endl;
    //Adjust match start
    new_match.shift_right(text_inspector_.vertex_left_siblings_length());
    //::std::cout << "After shift: " << new_match << ::std::endl;
    //cut the sequence to begin after the large_pattern_part_left_bound
    new_match.fit_into(first_lookup_begin_position_, last_lookup_begin_position_);
    //::std::cout << "After fit into [" << first_lookup_begin_position_ << ", " << last_lookup_begin_position_ << "]: " << new_match << ::std::endl;
    new_match.join_with(result);

    if (result &&
        (!new_match ||
          new_match.step() > pattern_.length() //any two consequent matches must have some common point
        )) {
      return result;
    }
    result = ::std::move(new_match);

    ++text_inspector_;
  }

  return result;
}

const FiniteArithmeticSequence& MatchingTable::matches(const Vertex& pattern,
                                                      const Vertex& text) {
  if (pattern.length() == 0) {
    return FiniteArithmeticSequence::Null;
  }

  if (pattern.length() > text.length()) {
    return FiniteArithmeticSequence::Null;
  }

  auto match_result_iterator = match_table_.find(std::make_pair(pattern, text));

  if (match_result_iterator != match_table_.end()) { //if already calculated
    return match_result_iterator->second;
  }

  FiniteArithmeticSequence match_result;

  if (pattern == text) { //If we are checking the same vertex
    match_result = FiniteArithmeticSequence(0, 1, 1);
  } else if (pattern.length() == 1) {//Trivial case
    Vertex pattern_vertex = pattern;
    while (pattern_vertex.height() > 1) {
      if (pattern_vertex.left_child() != Vertex::Null) {
        pattern_vertex = pattern_vertex.left_child();
      } else {
        pattern_vertex = pattern_vertex.right_child();
      }
    }

    Vertex text_letter_left_split = text; //the rightmost vertex in text.left_child()
    if (text.height() > 1) {
      text_letter_left_split = text.left_child();

      while (text_letter_left_split.height() > 1) {
        text_letter_left_split = text_letter_left_split.right_child() ? text_letter_left_split.right_child() : text_letter_left_split.left_child();
      }
    }


    Vertex text_letter_right_split = text.right_child(); //the leftmost vertex in text.right_child()
    while (text_letter_right_split.height() > 1) {
      text_letter_right_split = text_letter_right_split.left_child() ? text_letter_right_split.left_child() : text_letter_right_split.right_child();
    }

    if (pattern_vertex == text_letter_right_split && pattern_vertex == text_letter_left_split) {
      match_result = FiniteArithmeticSequence(text.left_child().length() - 1, 1, 2);
    } else if (pattern_vertex == text_letter_right_split) {
      match_result = FiniteArithmeticSequence(text.left_child().length(), 1, 1);
    } else if (pattern_vertex == text_letter_left_split) {
      match_result = FiniteArithmeticSequence(text.left_child().length() - 1, 1, 1);
    }
  } else  {//we have pattern.length > 1 => text.length > 1
    if (pattern.left_child().length() >= pattern.right_child().length()) {//Right child is smaller
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
    //::std::cout << "Not calculated match(";
    //PrintTo(pattern, &::std::cout);
    //::std::cout << ",";
    //PrintTo(text, &::std::cout);
    //::std::cout << ") is " << match_result << ::std::endl;

  }

  auto inserted_element = match_table_.insert(std::make_pair(std::make_pair(pattern, text), match_result));

  return inserted_element.first->second;
}

namespace internal {

FiniteArithmeticSequence nontrivial_match(
    const Vertex& large_pattern_part,
    const Vertex& small_pattern_part,
    bool small_pattern_is_after,
    const Vertex& text,
    MatchingTable* matching_table) {

  LongInteger large_pattern_part_left_bound = text.split_point() - large_pattern_part.length();

  if (small_pattern_is_after) {
    large_pattern_part_left_bound -= small_pattern_part.length();
  }

  PatternMatchesGenerator large_part_hunter(
      large_pattern_part, //pattern
      text, //text
      large_pattern_part_left_bound, //the leftmost letter of text to consider
      large_pattern_part.length() * 2 + small_pattern_part.length(),
      matching_table
  );

  FiniteArithmeticSequence result;

  while(large_part_hunter) {
    auto large_part_matches = large_part_hunter.next_match();

    //::std::cout << "Found large: " << large_part_matches << ::std::endl;

    if (large_part_matches) {
      FiniteArithmeticSequence small_part_candidates = large_part_matches;

      LongInteger seaside_candidates_length = 2 * small_pattern_part.length();
      if (small_pattern_is_after) {
        small_part_candidates.shift_right(large_pattern_part.length());
      } else {
        small_part_candidates.shift_right(-small_pattern_part.length());
      }
      //::std::cout << "Candidates: " << small_part_candidates << ::std::endl;
      LongInteger seaside_candidates_bound = (small_pattern_is_after ? small_part_candidates.last() - small_pattern_part.length(): small_part_candidates.first());

      PatternMatchesGenerator seaside_hunter(
          small_pattern_part,
          text,
          seaside_candidates_bound,
          2 * small_pattern_part.length(),
          matching_table
      );

      auto seaside_matches = seaside_hunter.next_match();
      //::std::cout << "Seaside matches: " << seaside_matches << ::std::endl;
      seaside_matches.intersect_with(small_part_candidates);
      //::std::cout << "Found seaside: " << seaside_matches << ::std::endl;

      const LongInteger& continental_candidate_start = small_pattern_is_after ? small_part_candidates.first() : small_part_candidates.last();
      PatternMatchesGenerator continental_hunter(
          small_pattern_part,
          text,
          continental_candidate_start,
          small_pattern_part.length(),
          matching_table
      );

      auto continental_matches = continental_hunter.next_match();

      FiniteArithmeticSequence& small_part_approved_candidates = seaside_matches;
      if (continental_matches) {
        //::std::cout << "Found continental: " << continental_matches << ::std::endl;
        continental_matches = small_part_candidates;
        if (small_pattern_is_after) {
          continental_matches.fit_into(continental_matches.first(), seaside_candidates_bound - 1);
        } else {
          continental_matches.fit_into(seaside_candidates_bound + small_pattern_part.length() + 1, continental_matches.last());
        }
        //::std::cout << "Multiplied continental: " << continental_matches << ::std::endl;


        small_part_approved_candidates.join_with(continental_matches);
      }


      if (small_pattern_is_after) {
        small_part_approved_candidates.shift_right(-large_pattern_part.length());
        small_part_approved_candidates.fit_into(0, text.length());
      }

      result.join_with(small_part_approved_candidates);
    }
  }

  return result;
}

} //namespace internal
} //namespace slp
} //namespace crag


