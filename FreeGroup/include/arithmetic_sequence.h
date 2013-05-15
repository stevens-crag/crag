/**
 * \file arithmetic_sequence.h
 * \brief Representation of a finite integer arithmetic sequence.
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
class FiniteArithmeticSequence {
  public:
    //! Default constructor
    FiniteArithmeticSequence()
      : first_(0)
      , step_(0)
      , last_(-1)
    { }

    //! Constructor taking first element, sequence step and elements count
    FiniteArithmeticSequence(LongInteger first, LongInteger step, LongInteger count)
      : first_(::std::move(first))
      , step_(::std::move(step))
      , last_(::std::move(((count -= 1) *= step_) += first_))
    {
      if (step_ <= 0 || last_ < first_) {
        first_ = step_ = 0;
        last_ = -1;
      } else if (last_ == first_) {
        step_ = 1;
      }
    }

    //All copy|move constructors can be defined by default

    //! Return first element
    const LongInteger& first() const {
      return first_;
    }

    //! Return last element (not after the last, but the last)
    const LongInteger& last() const {
      return last_;
    }

    //! Alias for last()
    const LongInteger& back() const {
      return last();
    }

    //! Return the difference between two sequence elements
    const LongInteger& step() const {
      return step_;
    }

    //! Return the number of elements. Note that it is calculated on each call.
    LongInteger count() const {
      return (last_ < first_) ? LongInteger(0) : ((last_ - first_) / step_ + 1);
    }

    //! Return the particular sequence element. No range check.
    LongInteger operator[](LongInteger index) const {
      return (index *= step_) + first_;
    }

    //! Check whether sequence contain some element
    bool contains(const LongInteger& position) const;

    //! Swap two sequences
    friend void swap(FiniteArithmeticSequence& first, FiniteArithmeticSequence& second) {
      using std::swap;

      swap(first.first_, second.first_);
      swap(first.step_, second.step_);
      swap(first.last_, second.last_);
    }

    //! False only if sequence is empty(), i.e. count() == 0
    explicit operator bool() const {
      return step_ != 0;
    }

    //! Check whether two sequences are the same
    bool operator==(const FiniteArithmeticSequence& other) const {
      return first_ == other.first_ && step_ == other.step_ && last_ == other.last_;
    }

    bool operator!=(const FiniteArithmeticSequence& other) const {
      return !(*this == other);
    }

    //! Shrink this so that for each element left_bound <= element <= right_bound
    FiniteArithmeticSequence& fit_into(const LongInteger& left_bound, const LongInteger& right_bound);

    //! Add \a shift to each element of sequence
    FiniteArithmeticSequence& shift_right(const LongInteger& shift) {
      if (*this) {
        first_ += shift;
        last_ += shift;
      }

      return *this;
    }

    //! Modify current sequence such that every element of the resulting sequence is element of this or other
    FiniteArithmeticSequence& join_with(const FiniteArithmeticSequence& other);
    //! Modify current sequence such that every element of the resulting sequence is element of this and other
    FiniteArithmeticSequence& intersect_with(const FiniteArithmeticSequence& other);

    //Use this method to get const reference to null sequence
    static const FiniteArithmeticSequence& NullSequence();
  private:
    LongInteger first_;
    LongInteger step_;
    LongInteger last_;
};

::std::ostream& operator<<(::std::ostream& out, const FiniteArithmeticSequence& sequence);

} //namespace crag


#endif /* CRAG_FREEGROUP_ARITHMETIC_SEQUENCE_H_ */
