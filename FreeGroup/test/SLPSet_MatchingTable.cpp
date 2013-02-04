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

TEST(MatchResultSequence, Swap) {
  SLPMatchingTable::MatchResultSequence first = {1, 2, 3};
  SLPMatchingTable::MatchResultSequence first_copy = {1, 2, 3};
  SLPMatchingTable::MatchResultSequence second = {4, 5, 6};
  SLPMatchingTable::MatchResultSequence second_copy = {4, 5, 6};

  std::swap(first, second);

  EXPECT_EQ(second_copy, first);
  EXPECT_EQ(first_copy, second);
}

TEST(JoinSequences, FirstZeroStep) {
  SLPMatchingTable::MatchResultSequence first = {1, 0, 1};
  SLPMatchingTable::MatchResultSequence second = {1, 2, 0};

  auto result = internal::join_sequences(first, second);

  EXPECT_EQ(result, second);
}

TEST(JoinSequences, SecondZeroStep) {
  SLPMatchingTable::MatchResultSequence first = {1, 1, 1};
  SLPMatchingTable::MatchResultSequence second = {1, 0, 2};

  auto result = internal::join_sequences(first, second);

  EXPECT_EQ(result, first);
}

TEST(JoinSequences, FirstZeroCount) {
  SLPMatchingTable::MatchResultSequence first = {1, 1, 0};
  SLPMatchingTable::MatchResultSequence second = {1, 2, 1};

  auto result = internal::join_sequences(first, second);

  EXPECT_EQ(result, second);
}

TEST(JoinSequences, SecondZeroCount) {
  SLPMatchingTable::MatchResultSequence first = {1, 1, 1};
  SLPMatchingTable::MatchResultSequence second = {1, 1, 0};

  auto result = internal::join_sequences(first, second);

  EXPECT_EQ(result, first);
}

TEST(JoinSequences, DifferentSteps) {
  SLPMatchingTable::MatchResultSequence first = {1, 1, 1};
  SLPMatchingTable::MatchResultSequence second = {1, 2, 1};

  auto result = internal::join_sequences(first, second);

  EXPECT_EQ(result, SLPMatchingTable::NO_MATCHES);
}

TEST(JoinSequences, NonCoherent) {
  SLPMatchingTable::MatchResultSequence first = {0, 2, 6};
  SLPMatchingTable::MatchResultSequence second = {1, 2, 10};
  SLPMatchingTable::MatchResultSequence expected = SLPMatchingTable::NO_MATCHES;

  auto result = internal::join_sequences(first, second);
  EXPECT_EQ(result, expected);
}

TEST(JoinSequences, NonIntersecting) {
  SLPMatchingTable::MatchResultSequence first = {0, 1, 1};
  SLPMatchingTable::MatchResultSequence second = {1, 1, 1};

  auto result = internal::join_sequences(first, second);

  EXPECT_EQ(result, SLPMatchingTable::NO_MATCHES);
}

TEST(JoinSequences, NonIntersectingSwapped) {
  SLPMatchingTable::MatchResultSequence first = {1, 1, 1};
  SLPMatchingTable::MatchResultSequence second = {0, 1, 1};

  auto result = internal::join_sequences(first, second);

  EXPECT_EQ(SLPMatchingTable::NO_MATCHES, result);
}

TEST(JoinSequences, InterectOneElement) {
  SLPMatchingTable::MatchResultSequence first = {0, 1, 2};
  SLPMatchingTable::MatchResultSequence second = {1, 1, 10};
  SLPMatchingTable::MatchResultSequence expected = {1, 1, 1};

  auto result = internal::join_sequences(first, second);
  EXPECT_EQ(result, expected);
}

TEST(JoinSequences, InterectOneElementSwapped) {
  SLPMatchingTable::MatchResultSequence first = {1, 1, 10};
  SLPMatchingTable::MatchResultSequence second = {0, 1, 2};
  SLPMatchingTable::MatchResultSequence expected = {1, 1, 1};

  auto result = internal::join_sequences(first, second);
  EXPECT_EQ(result, expected);
}

TEST(JoinSequences, InterectSeveralElements) {
  SLPMatchingTable::MatchResultSequence first = {0, 2, 6};
  SLPMatchingTable::MatchResultSequence second = {4, 2, 10};
  SLPMatchingTable::MatchResultSequence expected = {4, 2, 4};

  auto result = internal::join_sequences(first, second);
  EXPECT_EQ(result, expected);
}

TEST(JoinSequences, InterectSeveralElementsSwapped) {
  SLPMatchingTable::MatchResultSequence first = {4, 2, 10};
  SLPMatchingTable::MatchResultSequence second = {0, 2, 6};
  SLPMatchingTable::MatchResultSequence expected = {4, 2, 4};

  auto result = internal::join_sequences(first, second);
  EXPECT_EQ(result, expected);
}

TEST(JoinSequences, InterectSeveralElementsNonCoherentSwapped) {
  SLPMatchingTable::MatchResultSequence first = {0, 2, 6};
  SLPMatchingTable::MatchResultSequence second = {1, 2, 10};
  SLPMatchingTable::MatchResultSequence expected = SLPMatchingTable::NO_MATCHES;

  auto result = internal::join_sequences(first, second);
  EXPECT_EQ(result, expected);
}

TEST(JoinSequences, InterectSeveralElementsOneInsideAnother) {
  SLPMatchingTable::MatchResultSequence first = {0, 2, 10};
  SLPMatchingTable::MatchResultSequence second = {4, 2, 2};
  SLPMatchingTable::MatchResultSequence expected = {4, 2, 2};

  auto result = internal::join_sequences(first, second);
  EXPECT_EQ(result, expected);
}

TEST(JoinSequences, InterectSeveralElementsOneInsideAnotherSwapped) {
  SLPMatchingTable::MatchResultSequence first = {4, 2, 2};
  SLPMatchingTable::MatchResultSequence second = {0, 2, 10};
  SLPMatchingTable::MatchResultSequence expected = {4, 2, 2};

  auto result = internal::join_sequences(first, second);
  EXPECT_EQ(result, expected);
}

TEST(IntersectSequences, ZeroValues) {
  SLPMatchingTable::MatchResultSequence zero_step = {0, 0, 1};
  SLPMatchingTable::MatchResultSequence zero_count = {0, 1, 0};
  SLPMatchingTable::MatchResultSequence normal = {0, 1, 10};

  EXPECT_EQ(SLPMatchingTable::NO_MATCHES, internal::intersect_sequences(zero_step, normal));
  EXPECT_EQ(SLPMatchingTable::NO_MATCHES, internal::intersect_sequences(zero_count, normal));
}

TEST(IntersectSequences, SameSequences) {
  SLPMatchingTable::MatchResultSequence first = {10, 10, 10};
  SLPMatchingTable::MatchResultSequence second = {10, 10, 10};
  SLPMatchingTable::MatchResultSequence expected = {10, 10, 10};

  auto result = internal::intersect_sequences(first, second);
  EXPECT_EQ(expected, result);
}
TEST(IntersectSequences, InitialSupsequence) {
  SLPMatchingTable::MatchResultSequence first = {10, 10, 10};
  SLPMatchingTable::MatchResultSequence second = {10, 10, 2};
  SLPMatchingTable::MatchResultSequence expected = {10, 10, 2};

  auto result = internal::intersect_sequences(first, second);
  EXPECT_EQ(expected, result);
}
TEST(IntersectSequences, InitialSubsequence) {
  SLPMatchingTable::MatchResultSequence first = {10, 10, 2};
  SLPMatchingTable::MatchResultSequence second = {10, 10, 10};
  SLPMatchingTable::MatchResultSequence expected = {10, 10, 2};

  auto result = internal::intersect_sequences(first, second);
  EXPECT_EQ(expected, result);
}
TEST(IntersectSequences, TerminalSupsequence) {
  SLPMatchingTable::MatchResultSequence first = {10, 10, 10};
  SLPMatchingTable::MatchResultSequence second = {90, 10, 2};
  SLPMatchingTable::MatchResultSequence expected = {90, 10, 2};

  auto result = internal::intersect_sequences(first, second);
  EXPECT_EQ(expected, result);
}
TEST(IntersectSequences, TerminalSubsequence) {
  SLPMatchingTable::MatchResultSequence first = {90, 10, 2};
  SLPMatchingTable::MatchResultSequence second = {10, 10, 10};
  SLPMatchingTable::MatchResultSequence expected = {90, 10, 2};

  auto result = internal::intersect_sequences(first, second);
  EXPECT_EQ(expected, result);
}
TEST(IntersectSequences, InnerSupsequence) {
  SLPMatchingTable::MatchResultSequence first = {10, 10, 10};
  SLPMatchingTable::MatchResultSequence second = {50, 10, 2};
  SLPMatchingTable::MatchResultSequence expected = {50, 10, 2};

  auto result = internal::intersect_sequences(first, second);
  EXPECT_EQ(expected, result);
}
TEST(IntersectSequences, InnerSubsequence) {
  SLPMatchingTable::MatchResultSequence first = {50, 10, 2};
  SLPMatchingTable::MatchResultSequence second = {10, 10, 10};
  SLPMatchingTable::MatchResultSequence expected = {50, 10, 2};

  auto result = internal::intersect_sequences(first, second);
  EXPECT_EQ(expected, result);
}

TEST(IntersectSequences, InnerSupsequenceDifferentSteps) {
  SLPMatchingTable::MatchResultSequence first = {10, 5, 20};
  SLPMatchingTable::MatchResultSequence second = {50, 10, 2};
  SLPMatchingTable::MatchResultSequence expected = {50, 10, 2};

  auto result = internal::intersect_sequences(first, second);
  EXPECT_EQ(expected, result);
}
TEST(IntersectSequences, InnerSubsequenceDifferentSteps) {
  SLPMatchingTable::MatchResultSequence first = {50, 10, 2};
  SLPMatchingTable::MatchResultSequence second = {10, 5, 20};
  SLPMatchingTable::MatchResultSequence expected = {50, 10, 2};

  auto result = internal::intersect_sequences(first, second);
  EXPECT_EQ(expected, result);
}

TEST(IntersectSequences, NotCoherentStarts) {
  SLPMatchingTable::MatchResultSequence first = {10, 10, 10};
  SLPMatchingTable::MatchResultSequence second = {13, 15, 10};
  SLPMatchingTable::MatchResultSequence expected = SLPMatchingTable::NO_MATCHES;

  auto result = internal::intersect_sequences(first, second);
  EXPECT_EQ(expected, result);
}

TEST(IntersectSequences, OneAfterAnother) {
  SLPMatchingTable::MatchResultSequence first = {10, 10, 10};
  SLPMatchingTable::MatchResultSequence second = {0, 10, 2};
  SLPMatchingTable::MatchResultSequence expected = {10, 10, 1};

  auto result = internal::intersect_sequences(first, second);
  EXPECT_EQ(expected, result);
}

TEST(IntersectSequences, MiddleSingleIntersect) {
  SLPMatchingTable::MatchResultSequence first = {10, 10, 10};
  SLPMatchingTable::MatchResultSequence second = {35, 15, 2};
  SLPMatchingTable::MatchResultSequence expected = {50, 30, 1};

  auto result = internal::intersect_sequences(first, second);
  EXPECT_EQ(expected, result);
}

TEST(IntersectSequences, IntesectionOutOfBoundaries) {
  SLPMatchingTable::MatchResultSequence first = {55, 10, 10};
  SLPMatchingTable::MatchResultSequence second = {30, 15, 3};
  SLPMatchingTable::MatchResultSequence expected = SLPMatchingTable::NO_MATCHES;

  auto result = internal::intersect_sequences(first, second);
  EXPECT_EQ(expected, result);
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
      ASSERT_TRUE(current_first != current_second || current_first == current_intersection) <<
          "Common point of sequences " << first << " and " << second << " not in " << intersection;

      ASSERT_TRUE(intersection.count == 0 || (current_intersection >= current_first && current_intersection >= current_second)) <<
          "Extra point in intersection " << intersection << " of " << first << " and " << second;

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


TEST(FilledMatchingTable, Example) {
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

} //namespace
} //namespace crag
