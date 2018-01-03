/*
 * arithmetic_sequence.cpp
 *
 *  Created on: Feb 25, 2013
 *      Author: dpantele
 */

#include "gtest/gtest.h"

#include "arithmetic_sequence.h"

namespace crag{
namespace {

typedef FiniteArithmeticSequence Seq;

TEST(JoinArithmeticSequences, EmptySequences) {
  EXPECT_EQ(Seq(1, 2, 1), Seq(1, 0, 1).join_with(Seq(1, 2, 1)))
      << "If step of the first sequence is zero, return second.";

  EXPECT_EQ(Seq(1, 1, 1), Seq(1, 1, 1).join_with(Seq(1, 0, 2)))
      << "If step of the second sequence is zero, return first.";

  EXPECT_EQ(Seq(1, 2, 1), Seq(1, 1, 0).join_with(Seq(1, 2, 1)))
      << "If count of the first sequence is zero, return second.";

  EXPECT_EQ(Seq(1, 1, 1), Seq(1, 1, 1).join_with(Seq(1, 2, 0)))
      << "If count of the second sequence is zero, return first.";
}

TEST(JoinArithmeticSequences, StepCoherence) {
  EXPECT_EQ(Seq(), Seq(1, 1, 2).join_with(Seq(1, 2, 2)))
      << "If steps are different, the join is empty.";

  EXPECT_EQ(Seq(), Seq(0, 2, 6).join_with(Seq(1, 2, 10)))
      << "If the distance between first elements is not coherent with step, the join is empty.";
}

TEST(JoinArithmeticSequences, JoinTwoSequences) {
  EXPECT_EQ(Seq(), Seq(0, 2, 2).join_with(Seq(6, 2, 4)))
      << "First sequence ends far before the first element of the second.";

  EXPECT_EQ(Seq(0, 2, 7), Seq(0, 2, 3).join_with(Seq(6, 2, 4)))
      << "First sequence ends just before the first element of the second.";

  EXPECT_EQ(Seq(0, 2, 7), Seq(0, 2, 4).join_with(Seq(6, 2, 4)))
      << "First sequence ends at the first element of the second.";

  EXPECT_EQ(Seq(0, 2, 7), Seq(0, 2, 5).join_with(Seq(6, 2, 4)))
      << "First sequence starts before the second, ends after the first element of the second.";

  EXPECT_EQ(Seq(0, 2, 7), Seq(0, 2, 7).join_with(Seq(6, 2, 4)))
      << "First sequence starts before the second, ends at the last element of the second.";

  EXPECT_EQ(Seq(0, 2, 8), Seq(0, 2, 8).join_with(Seq(6, 2, 4)))
      << "First sequence starts before the second, ends after the last element of the second.";

  EXPECT_EQ(Seq(6, 2, 4), Seq(6, 2, 2).join_with(Seq(6, 2, 4)))
      << "First sequence starts with the second, ends after the first element of the second.";

  EXPECT_EQ(Seq(6, 2, 4), Seq(6, 2, 4).join_with(Seq(6, 2, 4)))
      << "Sequences are equal.";

  EXPECT_EQ(Seq(6, 2, 5), Seq(6, 2, 5).join_with(Seq(6, 2, 4)))
      << "First sequence starts with the second, ends after the last element of the second.";

  EXPECT_EQ(Seq(6, 2, 4), Seq(8, 2, 2).join_with(Seq(6, 2, 4)))
      << "First sequence starts after the second, ends before the last element of the second.";

  EXPECT_EQ(Seq(6, 2, 4), Seq(8, 2, 3).join_with(Seq(6, 2, 4)))
      << "First sequence starts after the second, ends at the last element of the second.";

  EXPECT_EQ(Seq(6, 2, 5), Seq(8, 2, 4).join_with(Seq(6, 2, 4)))
      << "First sequence starts after the second, ends after the last element of the second.";

  EXPECT_EQ(Seq(6, 2, 5), Seq(8, 2, 4).join_with(Seq(6, 2, 4)))
      << "First sequence starts after the second, ends after the last element of the second.";

  EXPECT_EQ(Seq(6, 2, 5), Seq(12, 2, 2).join_with(Seq(6, 2, 4)))
      << "First sequence starts at the last element of the second.";

  EXPECT_EQ(Seq(6, 2, 6), Seq(14, 2, 2).join_with(Seq(6, 2, 4)))
      << "First sequence starts just after the last element of the second.";

  EXPECT_EQ(Seq(), Seq(16, 2, 2).join_with(Seq(6, 2, 4)))
      << "First sequence starts far after the last element of the second.";

}

TEST(JoinArithmeticSequences, JoinTwoElements) {
  EXPECT_EQ(Seq(0, 10, 2), Seq(0, 3, 1).join_with(Seq(10, 7, 1))) <<
      "Two elements always can be joined. Here the first is before the second.";

  EXPECT_EQ(Seq(0, 10, 2), Seq(10, 3, 1).join_with(Seq(0, 7, 1))) <<
      "Two elements always can be joined. Here the second is before the first.";

  EXPECT_EQ(Seq(0, 1, 1), Seq(0, 3, 1).join_with(Seq(0, 7, 1))) <<
      "Two elements always can be joined. Here the second is equal to the first.";
}

TEST(JoinArithmeticSequences, JoinSequenceAndElement) {
  EXPECT_EQ(Seq(), Seq(0, 3, 1).join_with(Seq(4, 2, 3))) <<
      "Elements is before the sequence start, but is not in one step from it.";

  EXPECT_EQ(Seq(2, 2, 4), Seq(2, 3, 1).join_with(Seq(4, 2, 3))) <<
      "Elements is one step from the sequence start.";

  EXPECT_EQ(Seq(4, 2, 3), Seq(4, 3, 1).join_with(Seq(4, 2, 3))) <<
      "Elements is at the sequence start.";

  EXPECT_EQ(Seq(4, 2, 3), Seq(6, 3, 1).join_with(Seq(4, 2, 3))) <<
      "Elements is inside the sequence start.";

  EXPECT_EQ(Seq(), Seq(5, 3, 1).join_with(Seq(4, 2, 3))) <<
      "Elements is between the sequence first and last elements, but does not equal to any element.";

  EXPECT_EQ(Seq(4, 2, 3), Seq(8, 3, 1).join_with(Seq(4, 2, 3))) <<
      "Elements is the last sequence element.";

  EXPECT_EQ(Seq(4, 2, 4), Seq(10, 3, 1).join_with(Seq(4, 2, 3))) <<
      "Elements is one step from the last sequence element.";

  EXPECT_EQ(Seq(), Seq(12, 3, 1).join_with(Seq(4, 2, 3))) <<
      "Elements is after the last sequence element, but not in one step from it.";
}

TEST(IntersectArithmeticSequences, ZeroValues) {
  EXPECT_EQ(Seq(), Seq(0, 0, 1).intersect_with(Seq(0, 1, 10)));
  EXPECT_EQ(Seq(), Seq(0, 1, 0).intersect_with(Seq(0, 1, 10)));
}

TEST(IntersectArithmeticSequences, SubSequences) {
  EXPECT_EQ(Seq(10, 10, 10), Seq(0, 10, 11).intersect_with(Seq(10, 10, 10)));
  EXPECT_EQ(Seq(10, 10, 10), Seq(0, 10, 12).intersect_with(Seq(10, 10, 10)));
  EXPECT_EQ(Seq(10, 10, 10), Seq(10, 10, 10).intersect_with(Seq(10, 10, 11)));
  EXPECT_EQ(Seq(10, 10, 10), Seq(10, 10, 10).intersect_with(Seq(10, 10, 10)));
  EXPECT_EQ(Seq(10, 10, 10), Seq(10, 10, 11).intersect_with(Seq(10, 10, 10)));
  EXPECT_EQ(Seq(10, 10, 10), Seq(10, 10, 10).intersect_with(Seq(0, 10, 12)));
  EXPECT_EQ(Seq(10, 10, 10), Seq(10, 10, 10).intersect_with(Seq(0, 10, 11)));
}

TEST(IntersectArithmeticSequences, InnerSupsequenceDifferentSteps) {
  EXPECT_EQ(Seq(50, 10, 2), Seq(50, 10, 2).intersect_with(Seq(10, 5, 20)));
  EXPECT_EQ(Seq(50, 10, 2), Seq(10, 5, 20).intersect_with(Seq(50, 10, 2)));
}

TEST(IntersectArithmeticSequences, NotCoherentStarts) {
  EXPECT_EQ(Seq(), Seq(10, 10, 10).intersect_with(Seq(13, 15, 10)));
}

TEST(IntersectArithmeticSequences, ElementIntersection) {
  EXPECT_EQ(Seq(17, 1, 1), Seq(17, 13, 10).intersect_with(Seq(0, 17, 2)));
  EXPECT_EQ(Seq(50, 1, 1), Seq(10, 10, 10).intersect_with(Seq(35, 15, 2)));
}

TEST(IntersectArithmeticSequences, IntesectionOutOfBoundaries) {
  EXPECT_EQ(Seq(), Seq(55, 10, 10).intersect_with(Seq(30, 15, 3)));
}

TEST(IntersectArithmeticSequences, Example1) {
  EXPECT_EQ(Seq(0, 1, 1), Seq(0, 11, 2).intersect_with(Seq(0, 10, 2)));
}

TEST(IntersectArithmeticSequences, StressTest) {
  for (unsigned int test_code = 0; test_code < 01000000u /*8^6*/; ++test_code) {
    //We encode current test data using 3 bits for each of sequence parameters;
    Seq first(
        ((test_code >> 0) & 0x7u),       //start 0..7
        ((test_code >> 3) & 0x7u) + 10u, //step 10..17
        ((test_code >> 6) & 0x7u)        //count 0..7
    );

    Seq second(
        ((test_code >>  9) & 0x7u),       //start 0..7
        ((test_code >> 12) & 0x7u) + 10u, //step 10..17
        ((test_code >> 15) & 0x7u)        //count 0..7
    );

    Seq intersection = Seq(first).intersect_with(second);

    LongInteger current_first = first.first();
    LongInteger current_second = second.first();
    LongInteger current_intersection = intersection.first();

    LongInteger current_first_steps = 0;
    LongInteger current_second_steps = 0;
    LongInteger current_intersection_steps = 0;

    while(current_first <= first.last() && current_second <= second.last()) {
      ASSERT_TRUE(current_first != current_second || current_first >= current_intersection) <<
          "Common point of sequences " << first << " and " << second << " not in " << intersection;

      ASSERT_TRUE(current_first != current_second || current_first <= current_intersection) <<
          "Extra point " << current_intersection << " in intersection " << intersection << " of " << first << " and " << second;

      if (current_first == current_second) {

        ASSERT_LE(current_intersection, intersection.last()) <<
            "Intersection " << intersection << " of " << first << " and " << second << " is too short";

        current_intersection += intersection.step();
        ++current_intersection_steps;
      }

      if (current_first <= current_second) {
        current_first += first.step();
        ++current_first_steps;
      } else {
        current_second += second.step();
        ++current_second_steps;
      }
    }

    if (current_intersection_steps != 0) {
      ASSERT_EQ(intersection.last() + intersection.step(), current_intersection) <<
          "Intersection " << intersection << " of " << first << " and " << second << " is too long";
    }

  }
}

} //anonymous namespace
} //namespace crag
