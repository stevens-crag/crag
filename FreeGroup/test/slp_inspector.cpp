/*
 * SLPSet_inspector.cpp
 *
 *  Created on: Jan 22, 2013
 *      Author: dpantele
 */


#include <memory>
#include <functional>
#include <utility>
#include <vector>

#include "gtest/gtest.h"

#include "slp_inspector.h"

namespace crag {
namespace slp {
namespace {
  template <typename InspectorType>
  void check_inspector_path(InspectorType&& inspector, const ::std::vector<Vertex>& correct_path, const ::std::vector<LongInteger>& correct_length) {
    auto length = correct_length.begin();
    for (auto correct_vertex : correct_path) {
      ASSERT_FALSE(inspector.stopped());
      EXPECT_EQ(correct_vertex, inspector.vertex());
      EXPECT_EQ(*length, inspector.vertex_left_siblings_length());
      inspector.next();
      ++length;
    }
    ASSERT_TRUE(inspector.stopped());
  }


  TEST(PostorderInspector, EmptyInspector) {
    EXPECT_TRUE(PostorderInspector().stopped());
    EXPECT_TRUE(PostorderInspector().next().stopped());
  }

  TEST(PostorderInspector, SameChildren) {
    auto a = TerminalVertex(TerminalSymbol() + 1);
    auto a2 = NonterminalVertex(a, a);

    EXPECT_NO_FATAL_FAILURE(check_inspector_path(PostorderInspector(a2), ::std::vector<Vertex>({a, a, a2}), ::std::vector<LongInteger>({0, 1, 0})));
  }

  TEST(PostorderInspector, LeftTree) {
    auto a = TerminalVertex(TerminalSymbol() + 1);
    auto b = TerminalVertex(TerminalSymbol() + 2);
    auto c = TerminalVertex(TerminalSymbol() + 3);

    auto ab = NonterminalVertex(a, b);
    auto abc = NonterminalVertex(ab, c);

    EXPECT_NO_FATAL_FAILURE(check_inspector_path(PostorderInspector(abc), ::std::vector<Vertex>({a, b, ab, c, abc}), ::std::vector<LongInteger>({0, 1, 0, 2, 0})));
  }

  TEST(PostorderInspector, RightTree) {
    auto a = TerminalVertex(TerminalSymbol() + 1);
    auto b = TerminalVertex(TerminalSymbol() + 2);
    auto c = TerminalVertex(TerminalSymbol() + 3);

    auto bc = NonterminalVertex(b, c);
    auto abc = NonterminalVertex(a, bc);

    EXPECT_NO_FATAL_FAILURE(check_inspector_path(PostorderInspector(abc), ::std::vector<Vertex>({a, b, c, bc, abc}), ::std::vector<LongInteger>({0, 1, 2, 1, 0})));
  }

  TEST(PostorderInspector, CrossedTree) {
    auto a = TerminalVertex(TerminalSymbol() + 1);
    auto b = TerminalVertex(TerminalSymbol() + 2);

    auto ab = NonterminalVertex(a, b);
    auto ba_1 = NonterminalVertex(b, a.negate());

    auto abab_1 = NonterminalVertex(ab, ba_1.negate());

    EXPECT_NO_FATAL_FAILURE(check_inspector_path(
        PostorderInspector(abab_1),
        ::std::vector<Vertex>({a, b, ab, a, b.negate(), ba_1.negate(), abab_1}),
        ::std::vector<LongInteger>({0, 1, 0, 2, 3, 2, 0})
    ));
  }

  TEST(InorderInspector, EmptyInspector) {
    EXPECT_TRUE(InorderInspector().stopped());
    EXPECT_TRUE(InorderInspector().next().stopped());
  }

  TEST(InorderInspector, SameChildren) {
    auto a = TerminalVertex(TerminalSymbol() + 1);
    auto a2 = NonterminalVertex(a, a);

    EXPECT_NO_FATAL_FAILURE(check_inspector_path(InorderInspector(a2), ::std::vector<Vertex>({a, a2, a}), ::std::vector<LongInteger>({0, 0, 1})));
  }

  TEST(InorderInspector, LeftTree) {
    auto a = TerminalVertex(TerminalSymbol() + 1);
    auto b = TerminalVertex(TerminalSymbol() + 2);
    auto c = TerminalVertex(TerminalSymbol() + 3);

    auto ab = NonterminalVertex(a, b);
    auto abc = NonterminalVertex(ab, c);

    EXPECT_NO_FATAL_FAILURE(check_inspector_path(InorderInspector(abc), ::std::vector<Vertex>({a, ab, b, abc, c}), ::std::vector<LongInteger>({0, 0, 1, 0, 2})));
  }

  TEST(InorderInspector, RightTree) {
    auto a = TerminalVertex(TerminalSymbol() + 1);
    auto b = TerminalVertex(TerminalSymbol() + 2);
    auto c = TerminalVertex(TerminalSymbol() + 3);

    auto bc = NonterminalVertex(b, c);
    auto abc = NonterminalVertex(a, bc);

    EXPECT_NO_FATAL_FAILURE(check_inspector_path(InorderInspector(abc), ::std::vector<Vertex>({a, abc, b, bc, c}), ::std::vector<LongInteger>({0, 0, 1, 1, 2})));
  }

  TEST(InorderInspector, CrossedTree) {
    auto a = TerminalVertex(TerminalSymbol() + 1);
    auto b = TerminalVertex(TerminalSymbol() + 2);

    auto ab = NonterminalVertex(a, b);
    auto ba_1 = NonterminalVertex(b, a.negate());

    auto abab_1 = NonterminalVertex(ab, ba_1.negate());

    EXPECT_NO_FATAL_FAILURE(check_inspector_path(
        InorderInspector(abab_1),
        ::std::vector<Vertex>({a, ab, b, abab_1, a, ba_1.negate(), b.negate()}),
        ::std::vector<LongInteger>({0, 0, 1, 0, 2, 2, 3})
    ));
  }

  TEST(PreorderInspector, EmptyInspector) {
    EXPECT_TRUE(PreorderInspector().stopped());
    EXPECT_TRUE(PreorderInspector().next().stopped());
  }

  TEST(PreorderInspector, SameChildren) {
    auto a = TerminalVertex(TerminalSymbol() + 1);
    auto a2 = NonterminalVertex(a, a);

    EXPECT_NO_FATAL_FAILURE(check_inspector_path(PreorderInspector(a2), ::std::vector<Vertex>({a2, a, a}), ::std::vector<LongInteger>({0, 0, 1})));
  }

  TEST(PreorderInspector, LeftTree) {
    auto a = TerminalVertex(TerminalSymbol() + 1);
    auto b = TerminalVertex(TerminalSymbol() + 2);
    auto c = TerminalVertex(TerminalSymbol() + 3);

    auto ab = NonterminalVertex(a, b);
    auto abc = NonterminalVertex(ab, c);

    EXPECT_NO_FATAL_FAILURE(check_inspector_path(PreorderInspector(abc), ::std::vector<Vertex>({abc, ab, a, b, c}), ::std::vector<LongInteger>({0, 0, 0, 1, 2})));
  }

  TEST(PreorderInspector, RightTree) {
    auto a = TerminalVertex(TerminalSymbol() + 1);
    auto b = TerminalVertex(TerminalSymbol() + 2);
    auto c = TerminalVertex(TerminalSymbol() + 3);

    auto bc = NonterminalVertex(b, c);
    auto abc = NonterminalVertex(a, bc);

    EXPECT_NO_FATAL_FAILURE(check_inspector_path(PreorderInspector(abc), ::std::vector<Vertex>({abc, a, bc, b, c}), ::std::vector<LongInteger>({0, 0, 1, 1, 2})));
  }

  TEST(PreorderInspector, CrossedTree) {
    auto a = TerminalVertex(TerminalSymbol() + 1);
    auto b = TerminalVertex(TerminalSymbol() + 2);

    auto ab = NonterminalVertex(a, b);
    auto ba_1 = NonterminalVertex(b, a.negate());

    auto abab_1 = NonterminalVertex(ab, ba_1.negate());

    EXPECT_NO_FATAL_FAILURE(check_inspector_path(
        PreorderInspector(abab_1),
        ::std::vector<Vertex>({abab_1, ab, a, b, ba_1.negate(), a, b.negate()}),
        ::std::vector<LongInteger>({0, 0, 0, 1, 2, 2, 3})
    ));
  }
}
}
}


