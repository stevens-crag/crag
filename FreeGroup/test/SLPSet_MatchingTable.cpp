/*
 * SLPSet_MatchingTable.cpp
 *
 *  Created on: Jan 25, 2013
 *      Author: dpantele
 */

#include "gtest/gtest.h"
#include "SLPSet.h"

#include <memory>
#include <functional>
#include <utility>
#include <algorithm>
#include <cstdlib>
#include <ctime>


namespace crag {
::std::ostream& operator<<(::std::ostream& os, const ::crag::SLPVertex& vertex)
{
    return os << "vertex " << "(" << vertex.length().get_str() << ", " << vertex.height()
        << ", " << vertex.terminal_symbol() << ") - " << ::testing::PrintToString(vertex) << ::std::endl
        << "left " << "(" << vertex.left_child().length().get_str() << ", " << vertex.left_child().height()
        << ", " << vertex.left_child().terminal_symbol() << ") - " << ::testing::PrintToString(vertex.left_child()) << ::std::endl
        << "right " << "(" << vertex.right_child().length().get_str() << ", " << vertex.right_child().height()
        << ", " << vertex.right_child().terminal_symbol() << ") - " << ::testing::PrintToString(vertex.right_child()) << ::std::endl;
}

namespace {


class SLPMatchingInspectorTrivial : public ::testing::Test {
  protected:
    SLPVertex t;
    SLPVertex t2;

    SLPMatchingInspectorTrivial()
      : t(SLPVertex::terminal_vertex(1))
      , t2(SLPVertex::concatenate(t, t.negate()))
    { }
};

TEST_F(SLPMatchingInspectorTrivial, Empty) {
  internal::SLPMatchingInspector inspector(1, SLPVertex::Null, 0, 1);

  EXPECT_TRUE(inspector.inspection_ended())
    << "Nothing to visit anything";
}

TEST_F(SLPMatchingInspectorTrivial, OneVertexStraight) {
  internal::SLPMatchingInspector inspector(1, t, 0, 1);
  EXPECT_EQ(inspector.current_text(), t);
  inspector.go_next();
  EXPECT_TRUE(inspector.inspection_ended());
}

TEST_F(SLPMatchingInspectorTrivial, OneVertexLongPattern) {
  internal::SLPMatchingInspector inspector(2, t, 0, 2);
  EXPECT_TRUE(inspector.inspection_ended());
}

TEST_F(SLPMatchingInspectorTrivial, OneVertexShortCut) {
  internal::SLPMatchingInspector inspector(2, t, 0, 1);
  EXPECT_TRUE(inspector.inspection_ended());
}

TEST_F(SLPMatchingInspectorTrivial, OneVertexFarCut) {
  internal::SLPMatchingInspector inspector(1, t, 1, 1);
  EXPECT_TRUE(inspector.inspection_ended());
}

TEST_F(SLPMatchingInspectorTrivial, OneComposition) {
  internal::SLPMatchingInspector inspector(1, t2, 0, 10);
  EXPECT_EQ(inspector.current_text(), t);
  inspector.go_next();
  EXPECT_EQ(inspector.current_text(), t2);
  inspector.go_next();
  EXPECT_EQ(inspector.current_text(), t.negate());
  inspector.go_next();
  EXPECT_TRUE(inspector.inspection_ended());
}

TEST_F(SLPMatchingInspectorTrivial, OneCompositionLongPattern) {
  internal::SLPMatchingInspector inspector(2, t2, 0, 10);
  EXPECT_EQ(inspector.current_text(), t2);
  inspector.go_next();
  EXPECT_TRUE(inspector.inspection_ended());
}

TEST_F(SLPMatchingInspectorTrivial, OneCompositionLeftCut) {
  internal::SLPMatchingInspector inspector(1, t2, 0, 1);
  EXPECT_EQ(inspector.current_text(), t);
  inspector.go_next();
  EXPECT_EQ(inspector.current_text(), t2);
  inspector.go_next();
  EXPECT_TRUE(inspector.inspection_ended());
}

TEST_F(SLPMatchingInspectorTrivial, OneCompositionRightCut) {
  internal::SLPMatchingInspector inspector(1, t2, 1, 1);
  EXPECT_EQ(inspector.current_text(), t2);
  inspector.go_next();
  EXPECT_EQ(inspector.current_text(), t.negate());
  inspector.go_next();
  EXPECT_TRUE(inspector.inspection_ended());
}

typedef SLPMatchingTable::MatchResultSequence Seq;


TEST(MatchResultSequence, Swap) {
  Seq first = {1, 2, 3};
  Seq second = {4, 5, 6};

  std::swap(first, second);

  EXPECT_EQ(Seq({4, 5, 6}), first);
  EXPECT_EQ(Seq({1, 2, 3}), second);
}

using internal::join_sequences;
auto NO_MATCHES = SLPMatchingTable::NO_MATCHES;

TEST(JoinSequences, EmptySequnces) {
  EXPECT_EQ(Seq({1, 2, 1}), join_sequences({1, 0, 1}, {1, 2, 1}))
      << "If step of the first sequence is zero, return second.";

  EXPECT_EQ(Seq({1, 1, 1}), join_sequences({1, 1, 1}, {1, 0, 2}))
      << "If step of the second sequence is zero, return first.";

  EXPECT_EQ(Seq({1, 2, 1}), join_sequences({1, 1, 0}, {1, 2, 1}))
      << "If count of the first sequence is zero, return second.";

  EXPECT_EQ(Seq({1, 1, 1}), join_sequences({1, 1, 1}, {1, 2, 0}))
      << "If count of the second sequence is zero, return first.";
}

TEST(JoinSequences, StepCoherence) {
  EXPECT_EQ(NO_MATCHES, join_sequences({1, 1, 2}, {1, 2, 2}))
      << "If steps are different, the join is empty.";

  EXPECT_EQ(NO_MATCHES, join_sequences({0, 2, 6}, {1, 2, 10}))
      << "If the distance between first elements is not coherent with step, the join is empty.";
}

TEST(JoinSequences, JoinTwoSequnces) {
  EXPECT_EQ(NO_MATCHES, join_sequences({0, 2, 2}, {6, 2, 4}))
      << "First sequence ends far before the first element of the second.";

  EXPECT_EQ(Seq({0, 2, 7}), join_sequences({0, 2, 3}, {6, 2, 4}))
      << "First sequence ends just before the first element of the second.";

  EXPECT_EQ(Seq({0, 2, 7}), join_sequences({0, 2, 4}, {6, 2, 4}))
      << "First sequence ends at the first element of the second.";

  EXPECT_EQ(Seq({0, 2, 7}), join_sequences({0, 2, 5}, {6, 2, 4}))
      << "First sequence starts before the second, ends after the first element of the second.";

  EXPECT_EQ(Seq({0, 2, 7}), join_sequences({0, 2, 7}, {6, 2, 4}))
      << "First sequence starts before the second, ends at the last element of the second.";

  EXPECT_EQ(Seq({0, 2, 8}), join_sequences({0, 2, 8}, {6, 2, 4}))
      << "First sequence starts before the second, ends after the last element of the second.";

  EXPECT_EQ(Seq({6, 2, 4}), join_sequences({6, 2, 2}, {6, 2, 4}))
      << "First sequence starts with the second, ends after the first element of the second.";

  EXPECT_EQ(Seq({6, 2, 4}), join_sequences({6, 2, 4}, {6, 2, 4}))
      << "Sequences are equal.";

  EXPECT_EQ(Seq({6, 2, 5}), join_sequences({6, 2, 5}, {6, 2, 4}))
      << "First sequence starts with the second, ends after the last element of the second.";

  EXPECT_EQ(Seq({6, 2, 4}), join_sequences({8, 2, 2}, {6, 2, 4}))
      << "First sequence starts after the second, ends before the last element of the second.";

  EXPECT_EQ(Seq({6, 2, 4}), join_sequences({8, 2, 3}, {6, 2, 4}))
      << "First sequence starts after the second, ends at the last element of the second.";

  EXPECT_EQ(Seq({6, 2, 5}), join_sequences({8, 2, 4}, {6, 2, 4}))
      << "First sequence starts after the second, ends after the last element of the second.";

  EXPECT_EQ(Seq({6, 2, 5}), join_sequences({8, 2, 4}, {6, 2, 4}))
      << "First sequence starts after the second, ends after the last element of the second.";

  EXPECT_EQ(Seq({6, 2, 5}), join_sequences({12, 2, 2}, {6, 2, 4}))
      << "First sequence starts at the last element of the second.";

  EXPECT_EQ(Seq({6, 2, 6}), join_sequences({14, 2, 2}, {6, 2, 4}))
      << "First sequence starts just after the last element of the second.";

  EXPECT_EQ(NO_MATCHES, join_sequences({16, 2, 2}, {6, 2, 4}))
      << "First sequence starts far after the last element of the second.";

}

TEST(JoinSequences, JoinTwoElements) {
  EXPECT_EQ(Seq({0, 10, 2}), join_sequences({0, 3, 1}, {10, 7, 1})) <<
      "Two elements always can be joined. Here the first is before the second.";

  EXPECT_EQ(Seq({0, 10, 2}), join_sequences({10, 3, 1}, {0, 7, 1})) <<
      "Two elements always can be joined. Here the second is before the first.";

  EXPECT_EQ(Seq({0, 1, 1}), join_sequences({0, 3, 1}, {0, 7, 1})) <<
      "Two elements always can be joined. Here the second is equal to the first.";
}

TEST(JoinSequences, JoinSequenceAndElement) {
  EXPECT_EQ(NO_MATCHES, join_sequences({0, 3, 1}, {4, 2, 3})) <<
      "Elements is before the sequence start, but is not in one step from it.";

  EXPECT_EQ(Seq({2, 2, 4}), join_sequences({2, 3, 1}, {4, 2, 3})) <<
      "Elements is one step from the sequence start.";

  EXPECT_EQ(Seq({4, 2, 3}), join_sequences({4, 3, 1}, {4, 2, 3})) <<
      "Elements is at the sequence start.";

  EXPECT_EQ(Seq({4, 2, 3}), join_sequences({6, 3, 1}, {4, 2, 3})) <<
      "Elements is inside the sequence start.";

  EXPECT_EQ(NO_MATCHES, join_sequences({5, 3, 1}, {4, 2, 3})) <<
      "Elements is between the sequence first and last elements, but does not equal to any element.";

  EXPECT_EQ(Seq({4, 2, 3}), join_sequences({8, 3, 1}, {4, 2, 3})) <<
      "Elements is the last sequence element.";

  EXPECT_EQ(Seq({4, 2, 4}), join_sequences({10, 3, 1}, {4, 2, 3})) <<
      "Elements is one step from the last sequence element.";

  EXPECT_EQ(NO_MATCHES, join_sequences({12, 3, 1}, {4, 2, 3})) <<
      "Elements is after the last sequence element, but not in one step from it.";
}

using internal::intersect_sequences;

TEST(IntersectSequences, ZeroValues) {
  EXPECT_EQ(NO_MATCHES, intersect_sequences({0, 0, 1}, {0, 1, 10}));
  EXPECT_EQ(NO_MATCHES, intersect_sequences({0, 1, 0}, {0, 1, 10}));
}

TEST(IntersectSequences, SubSequnces) {
  EXPECT_EQ(Seq({10, 10, 10}), intersect_sequences({0, 10, 11}, {10, 10, 10}));
  EXPECT_EQ(Seq({10, 10, 10}), intersect_sequences({0, 10, 12}, {10, 10, 10}));
  EXPECT_EQ(Seq({10, 10, 10}), intersect_sequences({10, 10, 10}, {10, 10, 11}));
  EXPECT_EQ(Seq({10, 10, 10}), intersect_sequences({10, 10, 10}, {10, 10, 10}));
  EXPECT_EQ(Seq({10, 10, 10}), intersect_sequences({10, 10, 11}, {10, 10, 10}));
  EXPECT_EQ(Seq({10, 10, 10}), intersect_sequences({10, 10, 10}, {0, 10, 12}));
  EXPECT_EQ(Seq({10, 10, 10}), intersect_sequences({10, 10, 10}, {0, 10, 11}));
}

TEST(IntersectSequences, InnerSupsequenceDifferentSteps) {
  EXPECT_EQ(Seq({50, 10, 2}), intersect_sequences({50, 10, 2}, {10, 5, 20}));
  EXPECT_EQ(Seq({50, 10, 2}), intersect_sequences({10, 5, 20}, {50, 10, 2}));
}

TEST(IntersectSequences, NotCoherentStarts) {
  EXPECT_EQ(NO_MATCHES, intersect_sequences({10, 10, 10}, {13, 15, 10}));
}

TEST(IntersectSequences, ElementIntersection) {
  EXPECT_EQ(Seq({17, 1, 1}), intersect_sequences({17, 13, 10}, {0, 17, 2}));
  EXPECT_EQ(Seq({50, 1, 1}), intersect_sequences({10, 10, 10}, {35, 15, 2}));
}

TEST(IntersectSequences, IntesectionOutOfBoundaries) {
  EXPECT_EQ(NO_MATCHES, intersect_sequences({55, 10, 10}, {30, 15, 3}));
}

TEST(IntersectSequences, StressTest) {
  for (unsigned int test_code = 0; test_code < 01000000u /*8^6*/; ++test_code) {
    //We encode current test data using 3 bits for each of sequence parameters;
    SLPMatchingTable::MatchResultSequence first = {
        ((test_code >> 0) & 0x7u),       //start 0..7
        ((test_code >> 3) & 0x7u) + 10u, //step 10..17
        ((test_code >> 6) & 0x7u)        //count 0..7
    };

    SLPMatchingTable::MatchResultSequence second = {
        ((test_code >>  9) & 0x7u),       //start 0..7
        ((test_code >> 12) & 0x7u) + 10u, //step 10..17
        ((test_code >> 15) & 0x7u)        //count 0..7
    };

    SLPMatchingTable::MatchResultSequence intersection = internal::intersect_sequences(first, second);

    LongInteger current_first = first.start;
    LongInteger current_second = second.start;
    LongInteger current_intersection = intersection.start;

    LongInteger current_first_steps = 0;
    LongInteger current_second_steps = 0;
    LongInteger current_intersection_steps = 0;

    while(current_first_steps < first.count && current_second_steps < second.count) {
      ASSERT_TRUE(current_first != current_second || current_first >= current_intersection) <<
          "Common point of sequences " << first << " and " << second << " not in " << intersection;

      ASSERT_TRUE(current_first != current_second || current_first <= current_intersection) <<
          "Extra point " << current_intersection << " in intersection " << intersection << " of " << first << " and " << second;

      if (current_first == current_second) {
        ASSERT_LT(current_intersection_steps, intersection.count) <<
            "Intersection " << intersection << " of " << first << " and " << second << " is too short";

        current_intersection += intersection.step;
        ++current_intersection_steps;
      }

      if (current_first <= current_second) {
        current_first += first.step;
        ++current_first_steps;
      } else {
        current_second += second.step;
        ++current_second_steps;
      }
    }

    ASSERT_EQ(current_intersection_steps, intersection.count) <<
        "Intersection " << intersection << " of " << first << " and " << second << " is too long";

  }
}

TEST(SLPVertexHashtable, Test) {
  SLPVertex a = SLPVertex::terminal_vertex(1);
  SLPVertex acopy = a;
  SLPVertex a_ = a.negate();
  SLPVertex a2 = SLPVertex::concatenate(a, a);

  std::unordered_map<std::pair<SLPVertex, SLPVertex>, int> map;
  EXPECT_EQ(0, map.size());

  map.insert(std::make_pair(std::make_pair(a, a), 1));
  EXPECT_EQ(1, map.size());
  auto result = map.find(std::make_pair(a, a));
  ASSERT_NE(result, map.end());
  EXPECT_EQ(result->second, 1);

  map.insert(std::make_pair(std::make_pair(a, a), 1));
  EXPECT_EQ(1, map.size());
  result = map.find(std::make_pair(a, a));
  ASSERT_NE(result, map.end());
  EXPECT_EQ(result->second, 1);

  map.insert(std::make_pair(std::make_pair(a, a.negate()), 2));
  EXPECT_EQ(2, map.size());
  result = map.find(std::make_pair(a, a.negate()));
  ASSERT_NE(result, map.end());
  EXPECT_EQ(result->second, 2);

  map.insert(std::make_pair(std::make_pair(a.negate(), a), 3));
  EXPECT_EQ(3, map.size());
  result = map.find(std::make_pair(a.negate(), a));
  ASSERT_NE(result, map.end());
  EXPECT_EQ(result->second, 3);

  map.insert(std::make_pair(std::make_pair(a2, a), 4));
  EXPECT_EQ(4, map.size());
  result = map.find(std::make_pair(a2, a));
  ASSERT_NE(result, map.end());
  EXPECT_EQ(result->second, 4);

  map.insert(std::make_pair(std::make_pair(acopy, a), 5));
  EXPECT_EQ(4, map.size());
  result = map.find(std::make_pair(acopy, a));
  ASSERT_NE(result, map.end());
  EXPECT_EQ(result->second, 1);
}

std::string get_vertex_text_as_string(const SLPVertex& parent_vertex) {
  std::string presentation;

  for(const auto & current_vertex : SLPProducedWord(parent_vertex)) {
    presentation.push_back('a' - 1 + current_vertex.terminal_symbol());
  }

  return presentation;
}


class FilledMatchingTable : public SLPMatchingTable {
  public:
    FilledMatchingTable(const SLPVertex& pattern, const SLPVertex& text)
      : SLPMatchingTable()
    {
      SLPPostorderInspector text_inspector(text);

      while (!text_inspector.inspection_ended()) {
        std::string current_text;
        for(const auto & vertex : SLPProducedWord(text_inspector.current_vertex())) {
          current_text.push_back(vertex.is_negative() ? -vertex.terminal_symbol() : vertex.terminal_symbol());
        }

        SLPPostorderInspector pattern_inspector(pattern);
        while (!pattern_inspector.inspection_ended()) {
          MatchResultSequence result = NO_MATCHES;
          if (pattern_inspector.current_vertex().length() <= text_inspector.current_vertex().length() &&
              match_table_.find(std::make_pair(pattern_inspector.current_vertex(), text_inspector.current_vertex())) == match_table_.end()) {
            std::string current_pattern;
            for(const auto & vertex : SLPProducedWord(pattern_inspector.current_vertex())) {
              current_pattern.push_back(vertex.is_negative() ? -vertex.terminal_symbol() : vertex.terminal_symbol());
            }
            size_t current_match = text_inspector.current_vertex().split_point().get_ui() -
                pattern_inspector.current_vertex().length().get_ui();

            if (text_inspector.current_vertex().split_point() < pattern_inspector.current_vertex().length()) {
              current_match = 0;
            }

            size_t last_possible_match = text_inspector.current_vertex().split_point().get_ui();
//            std::cout << get_vertex_text_as_string(pattern_inspector.current_vertex()) << " in "
//                << get_vertex_text_as_string(text_inspector.current_vertex())
//                << " from " << current_match << " to " << last_possible_match << std::endl;

            current_match = current_text.find(current_pattern, current_match);
            while (current_match <= last_possible_match && current_match != std::string::npos) {
              if (result.count > 1) {
                result.count += 1;
              } else if (result.count == 1) {
                result.step = current_match - result.start;
                ++result.count;
              } else {
                result.start = current_match;
                result.step = 1;
                result.count = 1;
              }
              ++current_match;
              current_match = current_text.find(current_pattern, current_match);
            }
          }
          this->match_table_.insert(std::make_pair(std::make_pair(pattern_inspector.current_vertex(), text_inspector.current_vertex()), result));
          pattern_inspector.go_to_next_vertex();
        }
        text_inspector.go_to_next_vertex();
      }
    }
};


TEST(FilledMatchingTable, Example1) {
  SLPVertex a = SLPVertex::terminal_vertex(1);
  SLPVertex a2 = SLPVertex::concatenate(a, a);
  SLPVertex a4 = SLPVertex::concatenate(a2, a2);
  SLPVertex text = SLPVertex::concatenate(a4, a4);

  SLPVertex a3 = SLPVertex::concatenate(a, a2);
  SLPVertex pattern = SLPVertex::concatenate(a3, a2);

  FilledMatchingTable matching_table(pattern, text);

  EXPECT_EQ(FilledMatchingTable::MatchResultSequence({0, 1, 1}), matching_table.matches(a,a));
  EXPECT_EQ(FilledMatchingTable::NO_MATCHES, matching_table.matches(a2,a));
  EXPECT_EQ(FilledMatchingTable::NO_MATCHES, matching_table.matches(a3,a));
  EXPECT_EQ(FilledMatchingTable::NO_MATCHES, matching_table.matches(pattern,a));
  EXPECT_EQ(FilledMatchingTable::MatchResultSequence({0, 1, 2}), matching_table.matches(a,a2));
  EXPECT_EQ(FilledMatchingTable::MatchResultSequence({0, 1, 1}), matching_table.matches(a2,a2));
  EXPECT_EQ(FilledMatchingTable::NO_MATCHES, matching_table.matches(a3,a2));
  EXPECT_EQ(FilledMatchingTable::NO_MATCHES, matching_table.matches(pattern,a2));
  EXPECT_EQ(FilledMatchingTable::MatchResultSequence({1, 1, 2}), matching_table.matches(a,a4));
  EXPECT_EQ(FilledMatchingTable::MatchResultSequence({0, 1, 3}), matching_table.matches(a2,a4));
  EXPECT_EQ(FilledMatchingTable::MatchResultSequence({0, 1, 2}), matching_table.matches(a3,a4));
  EXPECT_EQ(FilledMatchingTable::NO_MATCHES, matching_table.matches(pattern,a4));
  EXPECT_EQ(FilledMatchingTable::MatchResultSequence({3, 1, 2}), matching_table.matches(a,text));
  EXPECT_EQ(FilledMatchingTable::MatchResultSequence({2, 1, 3}), matching_table.matches(a2,text));
  EXPECT_EQ(FilledMatchingTable::MatchResultSequence({1, 1, 4}), matching_table.matches(a3,text));
  EXPECT_EQ(FilledMatchingTable::MatchResultSequence({0, 1, 4}), matching_table.matches(pattern,text));
}

TEST(FilledMatchingTable, Example2) {
  SLPVertex a = SLPVertex::terminal_vertex(1);
  SLPVertex b = SLPVertex::terminal_vertex(2);

  SLPVertex ba = SLPVertex::concatenate(b, a);
  SLPVertex bab = SLPVertex::concatenate(ba, b);
  SLPVertex baba = SLPVertex::concatenate(bab, a);

  FilledMatchingTable table(ba, baba);
  EXPECT_EQ(Seq({2, 1, 1}), table.matches(ba, baba));
}

TEST(FilledMatchingTable, Example3) {
  SLPVertex a = SLPVertex::terminal_vertex(1);
  SLPVertex b = SLPVertex::terminal_vertex(2);

  SLPVertex ba = SLPVertex::concatenate(b, a);
  SLPVertex ab = SLPVertex::concatenate(a, b);
  SLPVertex baab = SLPVertex::concatenate(ba, ab);

  FilledMatchingTable table(baab, baab);
  EXPECT_EQ(Seq({0, 1, 1}), table.matches(baab, baab));
}


TEST(LocalSearch, Trivial) {
  SLPVertex a = SLPVertex::terminal_vertex(1);

  internal::SLPMatchingInspector inspector(a.length(), a, 0, 1);
  FilledMatchingTable matching_table(a, a);

  auto matches = internal::local_search(a, &inspector, &matching_table);
  EXPECT_EQ(SLPMatchingTable::MatchResultSequence({0, 1, 1}), matches);

  matches = internal::local_search(a, &inspector, &matching_table);
  EXPECT_EQ(SLPMatchingTable::NO_MATCHES, matches);
}

TEST(LocalSearch, LeftCut) {
  SLPVertex a = SLPVertex::terminal_vertex(1);
  SLPVertex a2 = SLPVertex::concatenate(a, a);

  internal::SLPMatchingInspector inspector(a.length(), a2, 0, 1);

  FilledMatchingTable matching_table(a, a2);

  auto matches = internal::local_search(a, &inspector, &matching_table);
  EXPECT_EQ(SLPMatchingTable::MatchResultSequence({0, 1, 1}), matches);

  matches = internal::local_search(a, &inspector, &matching_table);
  EXPECT_EQ(SLPMatchingTable::NO_MATCHES, matches);
}

TEST(LocalSearch, RightCut) {
  SLPVertex a = SLPVertex::terminal_vertex(1);
  SLPVertex a2 = SLPVertex::concatenate(a, a);

  internal::SLPMatchingInspector inspector(a.length(), a2, 1, 1);

  FilledMatchingTable matching_table(a, a2);

  auto matches = internal::local_search(a, &inspector, &matching_table);
  EXPECT_EQ(SLPMatchingTable::MatchResultSequence({1, 1, 1}), matches);

  matches = internal::local_search(a, &inspector, &matching_table);
  EXPECT_EQ(SLPMatchingTable::NO_MATCHES, matches);
}

TEST(LocalSearch, SimpleNontrivial) {
  SLPVertex a = SLPVertex::terminal_vertex(1);
  SLPVertex b = SLPVertex::terminal_vertex(2);
  SLPVertex ab = SLPVertex::concatenate(a, b);
  SLPVertex abab = SLPVertex::concatenate(ab, ab);

  internal::SLPMatchingInspector inspector(ab.length(), abab, 0, 4);

  FilledMatchingTable matching_table(ab, abab);

  auto matches = internal::local_search(ab, &inspector, &matching_table);
  EXPECT_EQ(Seq({0, 2, 2}), matches);

  matches = internal::local_search(ab, &inspector, &matching_table);
  EXPECT_EQ(NO_MATCHES, matches);
}

TEST(LocalSearch, SimpleNontrivialLeftCut) {
  SLPVertex a = SLPVertex::terminal_vertex(1);
  SLPVertex b = SLPVertex::terminal_vertex(2);
  SLPVertex ab = SLPVertex::concatenate(a, b);
  SLPVertex abab = SLPVertex::concatenate(ab, ab);

  internal::SLPMatchingInspector inspector(ab.length(), abab, 1, 4);

  FilledMatchingTable matching_table(ab, abab);

  auto matches = internal::local_search(ab, &inspector, &matching_table);
  EXPECT_EQ(Seq({2, 1, 1}), matches);

  matches = internal::local_search(ab, &inspector, &matching_table);
  EXPECT_EQ(NO_MATCHES, matches);
}

TEST(LocalSearch, SimpleNontrivialRightCut) {
  SLPVertex a = SLPVertex::terminal_vertex(1);
  SLPVertex b = SLPVertex::terminal_vertex(2);
  SLPVertex ab = SLPVertex::concatenate(a, b);
  SLPVertex abab = SLPVertex::concatenate(ab, ab);

  internal::SLPMatchingInspector inspector(ab.length(), abab, 0, 3);

  FilledMatchingTable matching_table(ab, abab);

  auto matches = internal::local_search(ab, &inspector, &matching_table);
  EXPECT_EQ(Seq({0, 1, 1}), matches);

  matches = internal::local_search(ab, &inspector, &matching_table);
  EXPECT_EQ(NO_MATCHES, matches);
}

TEST(LocalSearch, SimpleNontrivialSplitted) {
  SLPVertex a = SLPVertex::terminal_vertex(1);
  SLPVertex b = SLPVertex::terminal_vertex(2);
  SLPVertex ab = SLPVertex::concatenate(a, b);
  SLPVertex aba = SLPVertex::concatenate(ab, a);
  SLPVertex abaab = SLPVertex::concatenate(aba, ab);

  internal::SLPMatchingInspector ab_inspector(ab.length(), abaab, 0, 5);

  FilledMatchingTable ab_matching_table(ab, abaab);

  auto matches = internal::local_search(ab, &ab_inspector, &ab_matching_table);
  EXPECT_EQ(Seq({0, 1, 1}), matches);

  matches = internal::local_search(ab, &ab_inspector, &ab_matching_table);
  EXPECT_EQ(Seq({3, 1, 1}), matches);

  matches = internal::local_search(ab, &ab_inspector, &ab_matching_table);
  EXPECT_EQ(NO_MATCHES, matches);

  internal::SLPMatchingInspector a_inspector(a.length(), abaab, 0, 5);

  FilledMatchingTable a_matching_table(a, abaab);

  matches = internal::local_search(a, &a_inspector, &a_matching_table);
  EXPECT_EQ(Seq({0, 1, 1}), matches);

  matches = internal::local_search(a, &a_inspector, &a_matching_table);
  EXPECT_EQ(Seq({2, 1, 2}), matches);

  matches = internal::local_search(ab, &ab_inspector, &ab_matching_table);
  EXPECT_EQ(NO_MATCHES, matches);

}

SLPVertex get_random_slp_on_2_letters(unsigned int WORD_SIZE) {
  SLPVertex a = SLPVertex::terminal_vertex(1);
  SLPVertex b = SLPVertex::terminal_vertex(2);

  int random_word = rand() % (1 << WORD_SIZE);
  std::vector<unsigned int> random_word_split;
  for (unsigned int i = 1; i < WORD_SIZE; ++i) {
    random_word_split.push_back(i);
  }

  std::random_shuffle(random_word_split.begin(), random_word_split.end());

  std::vector<SLPVertex> word_presentation;
  for (unsigned int i = 0; i < WORD_SIZE; ++i) {
    word_presentation.push_back(((random_word & (1 << i)) ? b : a));
  }

  for (unsigned int split : random_word_split) {
    SLPVertex new_vertex = SLPVertex::concatenate(word_presentation[split - 1], word_presentation[split]);
    for (unsigned int i = split - new_vertex.left_child().length().get_ui(); i < split + new_vertex.right_child().length(); ++i) {
      word_presentation[i] = new_vertex;
    }
  }

  return word_presentation.front();
}

std::string get_vertex_height_in_inorder(const SLPVertex& vertex) {
  std::ostringstream split_string;
  internal::SLPMatchingInspector inspector(1, vertex, 0, vertex.length());
  while (!inspector.inspection_ended()) {
    split_string << inspector.current_text().height() << ",";

    inspector.go_next();
  }

  return split_string.str();
}

TEST(LocalSearch, RandomWord) {
  const unsigned int WORD_SIZE = 16;
  int REPEAT = 10000;

  srand(time(NULL));

  while (--REPEAT >= 0) {
    SLPVertex text = get_random_slp_on_2_letters(WORD_SIZE);

    int random_pattern_number = rand() % (2 * WORD_SIZE - 1);
    SLPPostorderInspector pattern_getter(text);
    int i = random_pattern_number;
    while (--i >= 0) {
      pattern_getter.go_to_next_vertex();
    }

    SLPVertex pattern = pattern_getter.current_vertex();

    unsigned int left_boundary = WORD_SIZE;
    unsigned int right_boundary = 0;
    while (left_boundary >= right_boundary) {
      left_boundary = rand() % WORD_SIZE;
      right_boundary = rand() % WORD_SIZE + 1;
    }

    unsigned int text_length = right_boundary - left_boundary;

    internal::SLPMatchingInspector text_inspector(pattern.length(), text, left_boundary, text_length);
    FilledMatchingTable matching_table(pattern, text);

    auto pattern_string = get_vertex_text_as_string(pattern);

    auto text_string = get_vertex_text_as_string(text);

    std::string split_string = get_vertex_height_in_inorder(text);

    unsigned int last_match_position = left_boundary;

    auto match = internal::local_search(pattern, &text_inspector, &matching_table);
    int last_approved_match = -1;
    while (match != NO_MATCHES) {
      int checked_matches = 0;
      unsigned int next_match_to_check = match.start.get_ui();
      while (checked_matches < match.count) {
        if (last_approved_match != next_match_to_check) {
          ASSERT_LE(next_match_to_check + pattern.length().get_ui(), right_boundary)
              << "Local search found match of " << pattern_string
              << " in " << text_string << "[" << left_boundary << ":" << right_boundary
              << "] out of boundaries. Split string: " << split_string
              << " pattern #" << random_pattern_number;
          last_match_position = text_string.find(pattern_string, last_match_position);
          ASSERT_NE(std::string::npos, last_match_position)
              << "Naive algorithm can't find any more matches of " << pattern_string
              << " in " << text_string << "[" << left_boundary << ":" << right_boundary
              << "], while local search found " << next_match_to_check
              << " Split string: " << split_string
              << " pattern #" << random_pattern_number;
          ASSERT_EQ(next_match_to_check, last_match_position)
              << "Naive algorithm can't find match " << next_match_to_check << " of " << pattern_string
              << " in " << text_string << "[" << left_boundary << ":" << right_boundary
              << "], found " << last_match_position << " instead. Split string: " << split_string
              << " pattern #" << random_pattern_number;
          ++last_match_position;
        }
        ++checked_matches;
        last_approved_match = next_match_to_check;
        next_match_to_check += match.step.get_ui();
      }
      match = internal::local_search(pattern, &text_inspector, &matching_table);
    }

    last_match_position = text_string.find(pattern_string, last_match_position);

    ASSERT_TRUE(last_match_position == std::string::npos || last_match_position + pattern.length().get_ui() > right_boundary)
      << "Naive algorithm found one more match of " << pattern_string
      << " in " << text_string << "[" << left_boundary << ":" << right_boundary
      << "] at position" << last_match_position
      << " Split string: " << split_string
      << " pattern #" << random_pattern_number;
  }
}

TEST(SLPMatchingTable, Trivial) {
  SLPVertex a = SLPVertex::terminal_vertex(1);
  SLPVertex a_copy = SLPVertex::terminal_vertex(1);
  SLPVertex a2 = SLPVertex::concatenate(a, a_copy);

  SLPMatchingTable matching_table;

  EXPECT_EQ(Seq({0, 1, 1}), matching_table.matches(a, a));
  EXPECT_EQ(Seq({0, 1, 1}), matching_table.matches(a, a_copy));
  EXPECT_EQ(Seq({0, 1, 2}), matching_table.matches(a, a2));
}

TEST(SLPMatchingTable, Example1) {
  SLPVertex a = SLPVertex::terminal_vertex(1);
  SLPVertex b = SLPVertex::terminal_vertex(2);
  SLPVertex ab = SLPVertex::concatenate(a, b);

  SLPMatchingTable matching_table;

  EXPECT_EQ(Seq({0, 1, 1}), matching_table.matches(a, ab));
  EXPECT_EQ(Seq({1, 1, 1}), matching_table.matches(b, ab));
  EXPECT_EQ(Seq({0, 1, 1}), matching_table.matches(ab, ab));
}

TEST(SLPMatchingTable, Example2) {
  SLPVertex a = SLPVertex::terminal_vertex(1);
  SLPVertex b = SLPVertex::terminal_vertex(2);
  SLPVertex ba = SLPVertex::concatenate(b, a);
  SLPVertex bb = SLPVertex::concatenate(b, b);
  SLPVertex bba = SLPVertex::concatenate(bb, a);

  SLPMatchingTable matching_table;

  EXPECT_EQ(NO_MATCHES, matching_table.matches(ba, bb));
  EXPECT_EQ(Seq({1, 1, 1}), matching_table.matches(ba, bba));
}

TEST(SLPMatchingTable, Example3) {
  SLPVertex a = SLPVertex::terminal_vertex(1);
  SLPVertex b = SLPVertex::terminal_vertex(2);

  SLPVertex ba = SLPVertex::concatenate(b, a);
  SLPVertex bab = SLPVertex::concatenate(ba, b);
  SLPVertex baba = SLPVertex::concatenate(bab, a);

  SLPMatchingTable table;
  EXPECT_EQ(Seq({2, 1, 1}), table.matches(ba, baba));
}

TEST(SLPMatchingTable, Example4) {
  SLPVertex a = SLPVertex::terminal_vertex(1);
  SLPVertex b = SLPVertex::terminal_vertex(2);

  SLPVertex bb = SLPVertex::concatenate(b, b);
  SLPVertex bbb = SLPVertex::concatenate(bb, b);
  SLPVertex abbb = SLPVertex::concatenate(a, bbb);

  SLPMatchingTable table;
  EXPECT_EQ(Seq({1, 1, 1}), table.matches(bb, abbb));
}

TEST(SLPMatchingTable, Example5) {
  SLPVertex a = SLPVertex::terminal_vertex(1);
  SLPVertex b = SLPVertex::terminal_vertex(2);

  SLPVertex aa = SLPVertex::concatenate(a, a);
  SLPVertex aab = SLPVertex::concatenate(aa, b);
  SLPVertex baab = SLPVertex::concatenate(b, aab);

  SLPMatchingTable table;
  EXPECT_EQ(Seq({0, 1, 1}), table.matches(aa, aa));
  EXPECT_EQ(Seq({0, 1, 1}), table.matches(aa, aab));
  EXPECT_EQ(Seq({1, 1, 1}), table.matches(aa, baab));
  EXPECT_EQ(Seq({1, 1, 1}), table.matches(aab, baab));
}

TEST(SLPMatchingTable, Example6) {
  SLPVertex a = SLPVertex::terminal_vertex(1);

  SLPVertex aa = SLPVertex::concatenate(a, a);
  SLPVertex aaa = SLPVertex::concatenate(aa, a);
  SLPVertex a5 = SLPVertex::concatenate(aa, aaa);
  SLPVertex a6 = SLPVertex::concatenate(a5, a);

  SLPMatchingTable table;
  EXPECT_EQ(Seq({0, 1, 1}), table.matches(aa, aa));
  EXPECT_EQ(Seq({0, 1, 2}), table.matches(aa, aaa));
  EXPECT_EQ(Seq({0, 1, 3}), table.matches(aa, a5));
  EXPECT_EQ(Seq({3, 1, 2}), table.matches(aa, a6));
  EXPECT_EQ(Seq({0, 1, 1}), table.matches(aaa, aaa));
  EXPECT_EQ(Seq({0, 1, 3}), table.matches(aaa, a5));
  EXPECT_EQ(Seq({2, 1, 2}), table.matches(aaa, a6));
  EXPECT_EQ(Seq({0, 1, 1}), table.matches(a5, a5));
  EXPECT_EQ(Seq({0, 1, 2}), table.matches(a5, a6));
}

TEST(SLPMatchingTable, InversedTest) {
  SLPVertex a = SLPVertex::terminal_vertex(1);
  SLPVertex b = SLPVertex::terminal_vertex(2);

  SLPVertex ab = SLPVertex::concatenate(a, b);
  SLPVertex a_1b = SLPVertex::concatenate(a.negate(), b);
  SLPVertex b_1a_1b = SLPVertex::concatenate(b.negate(), a_1b);

  SLPMatchingTable matching_table;

  EXPECT_EQ(Seq({1, 1, 1}), matching_table.matches(ab, b_1a_1b.negate()));
  EXPECT_EQ(Seq({0, 1, 1}), matching_table.matches(ab.negate(), b_1a_1b));
}

TEST(DISABLED_SLPMatchingTable, StressTest) {
  const unsigned int WORD_SIZE = 10;
  int repeat = 5000;

  srand(time(NULL));

  while (--repeat >= 0) {
    SLPVertex slp = get_random_slp_on_2_letters(WORD_SIZE);

    SLPMatchingTable real_matching_table;
    FilledMatchingTable computed_matching_table(slp, slp);

    SLPPostorderInspector text(slp);
    unsigned int text_vertex_number = 0;

    while (!text.inspection_ended()) {
      ++text_vertex_number;

      SLPPostorderInspector pattern(slp);
      unsigned int pattern_vertex_number = 0;
      while (!pattern.inspection_ended()) {
        ++pattern_vertex_number;
        ASSERT_EQ(computed_matching_table.matches(pattern.current_vertex(), text.current_vertex()), real_matching_table.matches(pattern.current_vertex(), text.current_vertex()))
          << "pattern #" << pattern_vertex_number << "(" << get_vertex_text_as_string(pattern.current_vertex()) << "," << get_vertex_height_in_inorder(pattern.current_vertex()) << ")"
          << ", text #" << text_vertex_number << "(" << get_vertex_text_as_string(text.current_vertex()) << "," << get_vertex_height_in_inorder(text.current_vertex()) << ")";
        pattern.go_to_next_vertex();
      }
      text.go_to_next_vertex();
    }
  }
}



} //namespace
} //namespace crag
