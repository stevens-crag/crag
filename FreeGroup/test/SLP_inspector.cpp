/*
 * SLPSet_inspector.cpp
 *
 *  Created on: Jan 22, 2013
 *      Author: dpantele
 */


#include "gtest/gtest.h"
#include "slp.h"

#include <memory>
#include <functional>
#include <utility>
#include <vector>

namespace crag {
namespace slp {
namespace {
  typedef TerminalVertexTemplate<char> TerminalVertex;

  template <typename InspectorType>
  void check_inspector_path(InspectorType&& inspector, const ::std::vector<Vertex>& correct_path) {
    for (auto correct_vertex : correct_path) {
      ASSERT_FALSE(inspector.stopped());
      EXPECT_EQ(correct_vertex, inspector.vertex());
      inspector.next();
    }
    ASSERT_TRUE(inspector.stopped());
  }


  TEST(PostorderInspector, EmptyInspector) {
    EXPECT_TRUE(PostorderInspector().stopped());
    EXPECT_TRUE(PostorderInspector().next().stopped());
  }

  TEST(PostorderInspector, SameChildren) {
    auto a = TerminalVertex('a');
    auto a2 = NonterminalVertex(a, a);

    EXPECT_NO_FATAL_FAILURE(check_inspector_path(PostorderInspector(a2), ::std::vector<Vertex>({a, a, a2})));
  }

  TEST(PostorderInspector, LeftTree) {
    auto a = TerminalVertex('a');
    auto b = TerminalVertex('b');
    auto c = TerminalVertex('c');

    auto ab = NonterminalVertex(a, b);
    auto abc = NonterminalVertex(ab, c);

    EXPECT_NO_FATAL_FAILURE(check_inspector_path(PostorderInspector(abc), ::std::vector<Vertex>({a, b, ab, c, abc})));
  }

  TEST(PostorderInspector, RightTree) {
    auto a = TerminalVertex('a');
    auto b = TerminalVertex('b');
    auto c = TerminalVertex('c');

    auto bc = NonterminalVertex(b, c);
    auto abc = NonterminalVertex(a, bc);

    EXPECT_NO_FATAL_FAILURE(check_inspector_path(PostorderInspector(abc), ::std::vector<Vertex>({a, b, c, bc, abc})));
  }

  TEST(PostorderInspector, CrossedTree) {
    auto a = TerminalVertex('a');
    auto b = TerminalVertex('b');

    auto ab = NonterminalVertex(a, b);
    auto ba_1 = NonterminalVertex(b, a.negate());

    auto abab_1 = NonterminalVertex(ab, ba_1.negate());

    EXPECT_NO_FATAL_FAILURE(check_inspector_path(PostorderInspector(abab_1), ::std::vector<Vertex>({a, b, ab, a, b.negate(), ba_1.negate(), abab_1})));
  }

  TEST(InorderInspector, EmptyInspector) {
    EXPECT_TRUE(InorderInspector().stopped());
    EXPECT_TRUE(InorderInspector().next().stopped());
  }

  TEST(InorderInspector, SameChildren) {
    auto a = TerminalVertex('a');
    auto a2 = NonterminalVertex(a, a);

    EXPECT_NO_FATAL_FAILURE(check_inspector_path(InorderInspector(a2), ::std::vector<Vertex>({a, a2, a})));
  }

  TEST(InorderInspector, LeftTree) {
    auto a = TerminalVertex('a');
    auto b = TerminalVertex('b');
    auto c = TerminalVertex('c');

    auto ab = NonterminalVertex(a, b);
    auto abc = NonterminalVertex(ab, c);

    EXPECT_NO_FATAL_FAILURE(check_inspector_path(InorderInspector(abc), ::std::vector<Vertex>({a, ab, b, abc, c})));
  }

  TEST(InorderInspector, RightTree) {
    auto a = TerminalVertex('a');
    auto b = TerminalVertex('b');
    auto c = TerminalVertex('c');

    auto bc = NonterminalVertex(b, c);
    auto abc = NonterminalVertex(a, bc);

    EXPECT_NO_FATAL_FAILURE(check_inspector_path(InorderInspector(abc), ::std::vector<Vertex>({a, abc, b, bc, c})));
  }

  TEST(InorderInspector, CrossedTree) {
    auto a = TerminalVertex('a');
    auto b = TerminalVertex('b');

    auto ab = NonterminalVertex(a, b);
    auto ba_1 = NonterminalVertex(b, a.negate());

    auto abab_1 = NonterminalVertex(ab, ba_1.negate());

    EXPECT_NO_FATAL_FAILURE(check_inspector_path(InorderInspector(abab_1), ::std::vector<Vertex>({a, ab, b, abab_1, a, ba_1.negate(), b.negate()})));
  }

  TEST(PreorderInspector, EmptyInspector) {
    EXPECT_TRUE(PreorderInspector().stopped());
    EXPECT_TRUE(PreorderInspector().next().stopped());
  }

  TEST(PreorderInspector, SameChildren) {
    auto a = TerminalVertex('a');
    auto a2 = NonterminalVertex(a, a);

    EXPECT_NO_FATAL_FAILURE(check_inspector_path(PreorderInspector(a2), ::std::vector<Vertex>({a2, a, a})));
  }

  TEST(PreorderInspector, LeftTree) {
    auto a = TerminalVertex('a');
    auto b = TerminalVertex('b');
    auto c = TerminalVertex('c');

    auto ab = NonterminalVertex(a, b);
    auto abc = NonterminalVertex(ab, c);

    EXPECT_NO_FATAL_FAILURE(check_inspector_path(PreorderInspector(abc), ::std::vector<Vertex>({abc, ab, a, b, c})));
  }

  TEST(PreorderInspector, RightTree) {
    auto a = TerminalVertex('a');
    auto b = TerminalVertex('b');
    auto c = TerminalVertex('c');

    auto bc = NonterminalVertex(b, c);
    auto abc = NonterminalVertex(a, bc);

    EXPECT_NO_FATAL_FAILURE(check_inspector_path(PreorderInspector(abc), ::std::vector<Vertex>({abc, a, bc, b, c})));
  }

  TEST(PreorderInspector, CrossedTree) {
    auto a = TerminalVertex('a');
    auto b = TerminalVertex('b');

    auto ab = NonterminalVertex(a, b);
    auto ba_1 = NonterminalVertex(b, a.negate());

    auto abab_1 = NonterminalVertex(ab, ba_1.negate());

    EXPECT_NO_FATAL_FAILURE(check_inspector_path(PreorderInspector(abab_1), ::std::vector<Vertex>({abab_1, ab, a, b, ba_1.negate(), a, b.negate()})));
  }
}
}
}


