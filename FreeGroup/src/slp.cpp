/*
 * slp.cpp
 *
 *  Created on: Feb 20, 2013
 *      Author: dpantele
 */


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

const inspector::InspectorTask inspector::InspectorTask::DO_NOTHING = inspector::InspectorTask();

}

::std::ostream& operator<<(::std::ostream& out, const FiniteAritmeticSequence& sequence) {
  return out << sequence.first() << " + " << sequence.step() << " * 0..." << sequence.count();
}

FiniteAritmeticSequence& FiniteAritmeticSequence::fit_into(const LongInteger& left_bound, const LongInteger& right_bound) {
  if (count_ == 0) {
    return *this;
  }

  //TODO: add static to speedup

  if (first_ < left_bound) {
    LongInteger count_reduce; //The number of steps which are outside of the boundary
    LongInteger correct_begin_offset; //Distance between boundary and first match
    LongInteger current_begin_offset = left_bound - first_;
    mpz_cdiv_qr(count_reduce.get_mpz_t(), correct_begin_offset.get_mpz_t(), current_begin_offset.get_mpz_t(), step_.get_mpz_t());
    count -= count_reduce;
    first += current_begin_offset;
    first -= correct_begin_offset;
  }

  LongInteger last_element = back();
  if (last_element > right_bound) {
    LongInteger count_reduce;
    LongInteger current_end_offset = last_element - right_bound;
    mpz_cdiv_q(count_reduce.get_mpz_t(), current_end_offset.get_mpz_t(), step_.get_mpz_t());
    count_ -= count_reduce;
  }
  return *this;
}


FiniteAritmeticSequence& FiniteAritmeticSequence::join_with(const FiniteAritmeticSequence& other) {
  if (!*this) {
    return *this = other;
  }

  if (!other) {
    return *this;
  }

  if (this->step_ != other.step_ && other.count_ != 1 && this->count_ != 1) {
    return *this = FiniteAritmeticSequence();
  }

  if (first_ > other.first_) {
    return *this = ::std::move(FiniteAritmeticSequence(other).join_with(*this));
  }

  //Let make this a bit messy
  static mpz_t distance_between_starts,
         steps_inside_starts_interval,
         residue;  //Should be 0
  //mpz_inits(distance_between_starts, steps_inside_starts_interval, residue, NULL);

  mpz_sub(distance_between_starts, other.first_.get_mpz_t(), this->first_.get_mpz_t());
  if (this->count_ == 1) {
    if (other.count_ == 1) {
      if (mpz_cmp_ui(distance_between_starts, 0) != 0) {
        mpz_set(this->step_.get_mpz_t(), distance_between_starts);
      }
    } else {
      this->step_ = other.step_;
    }
  }

  mpz_fdiv_qr(steps_inside_starts_interval, residue, distance_between_starts, step_.get_mpz_t());

  if (mpz_cmp_ui(residue, 0) != 0) { //starts are not coherent with step
    *this = FiniteAritmeticSequence();
  } else if (mpz_cmp(this->count_.get_mpz_t(), steps_inside_starts_interval) < 0) { //first sequence ends before second starts
    *this = FiniteAritmeticSequence();
  } else {
    mpz_add(steps_inside_starts_interval, steps_inside_starts_interval, other.count_.get_mpz_t());
    if (mpz_cmp(this->count_.get_mpz_t(), steps_inside_starts_interval) < 0) {//Second sequence ends after first
      mpz_swap(this->count_.get_mpz_t(), steps_inside_starts_interval);
    }
  }

  //mpz_clears(distance_between_starts, steps_inside_starts_interval, residue, NULL);
  return *this;
}

FiniteAritmeticSequence& FiniteAritmeticSequence::intersect_with(const FiniteAritmeticSequence& other) {
  if (!*this || !other) {
    return *this = FiniteAritmeticSequence();
  }

  if (other.first_ > this->first_) { //Ensure that this starts after other
    return *this = ::std::move(FiniteAritmeticSequence(other).intersect_with(*this));
  }
  //from this point we assign index 1 to the other sequence and index 2 to this sequence, so first_1 <= first_2

  static mpz_t temp, //static here gives 5% performance
        step_2_inverse,
        steps_gcd,
        intersection_end;
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
    *this = FiniteAritmeticSequence();
  } else {
    mpz_mul(step_2_inverse, step_2_inverse, temp); // step_2_inverse = (first_1 - first_2) * step_2_inverse

    //now calculate the right boundary for the resulting sequence
    //first, take first_2 + step_2 * count_2 as an assumption (it starts after, maybe it also ends after)
    mpz_set(intersection_end, this->first_.get_mpz_t());
    mpz_addmul(intersection_end, this->step_.get_mpz_t(), this->count_.get_mpz_t());

    //temp is the end of this sequence #1 = first_1 + step_1 * count_1
    mpz_set(temp, other.first_.get_mpz_t());
    mpz_addmul(temp, other.step_.get_mpz_t(), other.count_.get_mpz_t());

    if (mpz_cmp(intersection_end, temp) > 0) { //take minimum of the ends of two sequences as the right boundary
      mpz_swap(intersection_end, temp);
    }

    //we store step_2 / gcd in this->step_
    mpz_divexact(this->step_.get_mpz_t(), this->step_.get_mpz_t(), steps_gcd);

    //now we have to calculate the start of resulting sequence
    mpz_fdiv_r(temp, step_2_inverse, other.step_.get_mpz_t()); //temp = k * gcd = step_2_inverse % step_1
    mpz_addmul(this->first_.get_mpz_t(), temp, this->step_.get_mpz_t()); //first_result = first_2 + (k * gcd) * (step_2 / gcd)

    //now calculate the result step as (step_2 / gcd) * step_1
    this->step_ *= other.step_;

    //now we just calculate the number of steps inside result just as a floor of (intersection_end - first_result)/step_result
    mpz_sub(temp, intersection_end, this->first_.get_mpz_t());
    mpz_cdiv_q(this->count_.get_mpz_t(), temp, this->step_.get_mpz_t());

    if (this->count_ <= 0) {
      *this = FiniteAritmeticSequence();
    } else if (this->count_ == 1) {
      this->step_ = 1;
    }
  }

  //mpz_clears(temp, step_2_inverse, steps_gcd, intersection_end, NULL);

  return *this;
}

const FiniteAritmeticSequence FiniteAritmeticSequence::Null = FiniteAritmeticSequence();

namespace slp {

const FiniteAritmeticSequence& PatternMatchesGenerator::next_match() {
  current_match_ = FiniteAritmeticSequence::Null;
  FiniteAritmeticSequence new_match;

  while (!text_inspector_.stopped()) {
    new_match = matching_table_->matches(pattern_, text_inspector_.vertex());
    //Adjust match start
    new_match.shift_right(text_inspector_.vertex_left_siblings_length());
    //cut the sequence to begin after the large_pattern_part_left_bound
    new_match.fit_into(first_lookup_begin_position_, last_lookup_begin_position_);
    new_match.join_with(current_match_);

    if (current_match_ &&
        (!new_match ||
          new_match.step > pattern_.length() //any two consequent matches must have some common point
        )) {
      return current_match_;
    }
    current_match_ = ::std::move(new_match);

    ++text_inspector_;
  }

  return current_match_;
}

} //namespace slp

} //namespace crag


