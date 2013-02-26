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

class FiniteAritmeticSequence {
  public:
    FiniteAritmeticSequence()
      : first_(0)
      , step_(0)
      , count_(0)
    { }

    FiniteAritmeticSequence(LongInteger first, LongInteger step, LongInteger count)
      : first_(::std::move(first))
      , step_(::std::move(step))
      , count_(::std::move(count))
    {
      if (count_ < 0 || step_ <= 0) {
        count_ = 0;
      }

      if (count_ == 0) {
        first_ = 0;
        step_ = 0;
      } else if (count_ == 1) {
        step_ = 1;
      }
    }

    //All copy|move constructors can be defined by default

    const LongInteger& first() const {
      return first_;
    }

    LongInteger back() const {
      return (*this)[count_ - 1];
    }

    const LongInteger& size() const {
      return count_;
    }

    const LongInteger& count() const {
      return size();
    }

    LongInteger operator[](const LongInteger& index) const {
      return first_ + index * step_;
    }

    const LongInteger& step() const {
      return step_;
    }

    friend void swap(FiniteAritmeticSequence& first, FiniteAritmeticSequence& second) {
      using std::swap;

      swap(first.first_, second.first_);
      swap(first.step_, second.step_);
      swap(first.count_, second.count_);
    }

    explicit operator bool() const {
      return count_ != 0;
    }

    bool operator==(const FiniteAritmeticSequence& other) const {
      return first_ == other.first_ && step_ == other.step_ && count_ == other.count_;
    }

    bool operator!=(const FiniteAritmeticSequence& other) const {
      return !(*this == other);
    }

    FiniteAritmeticSequence& join_with(const FiniteAritmeticSequence& other);
    FiniteAritmeticSequence& intersect_with(const FiniteAritmeticSequence& other);
  private:
    LongInteger first_;
    LongInteger step_;
    LongInteger count_;
};

::std::ostream& operator<<(::std::ostream& out, const FiniteAritmeticSequence& sequence);

} //namespace crag


#endif /* CRAG_FREEGROUP_ARITHMETIC_SEQUENCE_H_ */
