/*
 * SLPSet_inspector.cpp
 *
 *  Created on: Jan 22, 2013
 *      Author: dpantele
 */


#include "gtest/gtest.h"
#include "slp_inspector.h"
#include "slp_pattern_matching.h"

#include <memory>
#include <functional>
#include <utility>
#include <vector>
#include <tuple>

namespace crag {
namespace slp {
namespace {
  typedef TerminalVertexTemplate<char> TerminalVertex;

  typedef ::std::tuple<int, int, int, int, ::std::vector<int>> BoundedTaskAcceptorTestParam;

  class BoundedTaskAcceptorTest : public ::testing::TestWithParam<BoundedTaskAcceptorTestParam>
  {  };

  class PatternMatchesGeneratorAccess : public PatternMatchesGenerator {
    public:
      PatternMatchesGeneratorAccess(const Vertex& pattern, const Vertex& text, LongInteger lookup_from, LongInteger lookup_length)
        : PatternMatchesGenerator(pattern, text, ::std::move(lookup_from), ::std::move(lookup_length))
      { }

      InorderInspector& get_inspector() {
        return text_inspector_;
      }
  };

  TEST_P(BoundedTaskAcceptorTest, CheckOrder) {
    int pattern_code;
    int text_code;
    int lookup_from_int;
    int lookup_length_int;
    ::std::vector<int> correct_path;

    ::std::tie(pattern_code, text_code, lookup_from_int, lookup_length_int, correct_path) = GetParam();

    TerminalVertex t('t');
    NonterminalVertex t2(t, t);

    Vertex pattern = Vertex::Null, text = Vertex::Null;
    switch(pattern_code) {
    case 1:
      pattern = t;
      break;
    case 2:
      pattern = t2;
      break;
    }

    switch(text_code) {
    case 1:
      text = t;
      break;
    case 2:
      text = t2;
      break;
    }

    //Must update this code if you change constructors in PatternMatchesGenerator
    LongInteger lookup_from = lookup_from_int;
    LongInteger lookup_length = lookup_length_int;
    LongInteger first_lookup_begin_position_ = lookup_from;
    LongInteger first_lookup_end_position_(lookup_from + pattern.length());
    LongInteger last_lookup_begin_position_(((lookup_from += lookup_length) > text.length()) ? (text.length() - pattern.length()) : (lookup_from - pattern.length()));
    InorderInspector text_inspector_(text, ::std::unique_ptr<inspector::BoundedTaskAcceptor>(new inspector::BoundedTaskAcceptor(
      first_lookup_end_position_, last_lookup_begin_position_
    )));

    InorderInspector& text_inspector = text_inspector_;

    for (auto correct_vertex : correct_path) {
      ASSERT_FALSE(text_inspector.stopped());
      if (correct_vertex == 1) {
        EXPECT_EQ(t, text_inspector.vertex());
      } else {
        EXPECT_EQ(t2, text_inspector.vertex());
      }

      text_inspector.next();
    }
    ASSERT_TRUE(text_inspector.stopped());
  }

  INSTANTIATE_TEST_CASE_P(
      BoundedTaskAcceptor,
      BoundedTaskAcceptorTest,
      ::testing::Values(
        ::std::make_tuple(1, 0, 0, 1, ::std::vector<int>()),
        ::std::make_tuple(1, 1, 0, 1, ::std::vector<int>({1})),
        ::std::make_tuple(2, 1, 0, 2, ::std::vector<int>({})),
        ::std::make_tuple(2, 1, 0, 1, ::std::vector<int>({})),
        ::std::make_tuple(1, 1, 1, 1, ::std::vector<int>({})),
        ::std::make_tuple(1, 2, 0, 10, ::std::vector<int>({1, 2, 1})),
        ::std::make_tuple(2, 2, 0, 10, ::std::vector<int>({2})),
        ::std::make_tuple(1, 2, 0, 1, ::std::vector<int>({1, 2})),
        ::std::make_tuple(1, 2, 1, 1, ::std::vector<int>({2, 1}))
      ));

}
}
}


