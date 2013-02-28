/*
 * arithmetic_sequence.h
 *
 *  Created on: Feb 24, 2013
 *      Author: dpantele
 */

#pragma once
#ifndef CRAG_FREEGROUP_ARITHMETIC_SEQUENCE_H_
#define CRAG_FREEGROUP_ARITHMETIC_SEQUENCE_H_

#include <iostream>
#include <utility>

#include <gmpxx.h>
typedef mpz_class LongInteger; //Note that actually we are tied to GMP, because we use its internal a lot

namespace crag {

//! Finite arithmetic sequence
class FiniteAritmeticSequence {
  public:
    FiniteAritmeticSequence()
      : first_(0)
      , step_(0)
      , last_(-1)
    { }

    FiniteAritmeticSequence(LongInteger first, LongInteger step, LongInteger count)
      : first_(::std::move(first))
      , step_(::std::move(step))
      , last_(::std::move(((count -= 1) *= step) += first))
    {
      if (step_ <= 0 || last_ < first_) {
        first_ = step_ = 0;
        last_ = -1;
      } else if (last_ == first_) {
        step_ = 1;
      }
    }

    //All copy|move constructors can be defined by default

    const LongInteger& first() const {
      return first_;
    }

    const LongInteger& last() const {
      return last_;
    }

    const LongInteger& back() const {
      return last();
    }

    const LongInteger& step() const {
      return step_;
    }

    LongInteger count() const {
      return (last_ < first_) ? LongInteger(0) : ((last_ - first_) / step_);
    }

    LongInteger operator[](LongInteger index) const {
      return (index *= step_) + first_;
    }

    friend void swap(FiniteAritmeticSequence& first, FiniteAritmeticSequence& second) {
      using std::swap;

      swap(first.first_, second.first_);
      swap(first.step_, second.step_);
      swap(first.last_, second.last_);
    }

    explicit operator bool() const {
      return step_ != 0;
    }

    bool operator==(const FiniteAritmeticSequence& other) const {
      return first_ == other.first_ && step_ == other.step_ && last_ == other.last_;
    }

    bool operator!=(const FiniteAritmeticSequence& other) const {
      return !(*this == other);
    }

    FiniteAritmeticSequence& fit_into(const LongInteger& left_bound, const LongInteger& right_bound);
    FiniteAritmeticSequence& shift_right(const LongInteger& length) {
      if (*this) {
        first_ += length;
        last_ += length;
      }

      return *this;
    }

    FiniteAritmeticSequence& join_with(const FiniteAritmeticSequence& other);
    FiniteAritmeticSequence& intersect_with(const FiniteAritmeticSequence& other);

    const static FiniteAritmeticSequence Null;
  private:
    LongInteger first_;
    LongInteger step_;
    LongInteger last_;
};

::std::ostream& operator<<(::std::ostream& out, const FiniteAritmeticSequence& sequence);

} //namespace crag


#endif /* CRAG_FREEGROUP_ARITHMETIC_SEQUENCE_H_ */
