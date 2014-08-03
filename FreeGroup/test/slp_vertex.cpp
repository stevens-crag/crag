/*
 * SLPSet_common.cpp
 *
 *  Created on: Jan 21, 2013
 *      Author: dpantele
 */

#include <memory>
#include <functional>
#include <utility>

#include "gtest/gtest.h"
#include "slp_vertex.h"

namespace crag {
namespace slp {
namespace {

CONSTEXPR_OR_CONST TerminalSymbol terminal_a = TerminalSymbol{} + 100;
CONSTEXPR_OR_CONST TerminalSymbol terminal_b = TerminalSymbol{} + 200;
CONSTEXPR_OR_CONST TerminalSymbol terminal_c = TerminalSymbol{} + 300;

class VertexTest : public ::testing::Test {
  protected:
    VertexTest()
      : a(TerminalVertex(terminal_a))
      , b(TerminalVertex(terminal_b))
      , c(TerminalVertex(terminal_c))
    { }

    Vertex a;
    Vertex b;
    Vertex c;
    ::std::hash<Vertex> hash;
};

TEST_F(VertexTest, NullVertex) {
  EXPECT_EQ(Vertex(), Vertex().negate());
  EXPECT_EQ(Vertex(), Vertex(Vertex()));

  EXPECT_EQ(0, Vertex().length());
  EXPECT_EQ(0, Vertex().height());
  EXPECT_EQ(0, Vertex().split_point());
  EXPECT_EQ(Vertex(), Vertex().left_child());
  EXPECT_EQ(Vertex(), Vertex().right_child());

  EXPECT_EQ(std::hash<Vertex::VertexSignedId>()(0), hash(Vertex()));
  EXPECT_EQ(std::hash<Vertex::VertexSignedId>()(0), hash(Vertex().negate()));
}

TEST_F(VertexTest, TerminalVertex) {
  EXPECT_EQ(1, a.length());
  EXPECT_EQ(1, a.negate().length());
  EXPECT_EQ(1, a.height());
  EXPECT_EQ(1, a.negate().height());
  EXPECT_EQ(0, a.split_point());
  EXPECT_EQ(0, a.negate().split_point());

  EXPECT_EQ(0, Vertex().split_point());
  EXPECT_EQ(Vertex(), a.left_child());
  EXPECT_EQ(Vertex(), a.right_child());

  EXPECT_NE(a, b);
  EXPECT_NE(a, Vertex());
  EXPECT_NE(a, a.negate());
  EXPECT_EQ(a, a);
  EXPECT_EQ(a, a.negate().negate());
  EXPECT_EQ(terminal_a, TerminalVertex(a));
  EXPECT_EQ(-terminal_a, TerminalVertex(a.negate()).terminal_symbol());

  TerminalVertex typed_a(a);
  EXPECT_EQ(terminal_a, typed_a);
  EXPECT_EQ(-terminal_a, TerminalVertex(typed_a.negate()));

  EXPECT_EQ(hash(a), hash(a.negate().negate()));
  EXPECT_NE(hash(a), hash(a.negate()));
  EXPECT_NE(hash(a), hash(b));
  EXPECT_NE(hash(a), hash(Vertex()));
}

TEST_F(VertexTest, NonterminalVertex) {
  auto ab_1 = NonterminalVertex(a, b.negate());
  EXPECT_EQ(2, ab_1.length());
  EXPECT_EQ(2, ab_1.height());
  EXPECT_EQ(1, ab_1.split_point());
  EXPECT_EQ(2, ab_1.negate().length());
  EXPECT_EQ(2, ab_1.negate().height());
  EXPECT_EQ(1, ab_1.negate().split_point());
  EXPECT_EQ(a, ab_1.left_child());
  EXPECT_NE(b, ab_1.right_child());
  EXPECT_EQ(b.negate(), ab_1.right_child());

  auto ba_1c_1  = NonterminalVertex(ab_1.negate(), c.negate()); //ba^-1c^-11
  EXPECT_EQ(3, ba_1c_1.length());
  EXPECT_EQ(3, ba_1c_1.height());
  EXPECT_EQ(ab_1, ba_1c_1.left_child().negate());
  EXPECT_EQ(b, ba_1c_1.left_child().left_child());
  EXPECT_EQ(a.negate(), ba_1c_1.left_child().right_child());
  EXPECT_EQ(c.negate(), ba_1c_1.right_child());

  EXPECT_EQ(0, TerminalVertex(ab_1).terminal_symbol());

  EXPECT_EQ(hash(ab_1), hash(ab_1));
  EXPECT_EQ(hash(ab_1), hash(ba_1c_1.left_child().negate()));
  EXPECT_NE(hash(ab_1), hash(ab_1.negate()));
  EXPECT_NE(hash(ab_1), hash(ba_1c_1));
  EXPECT_NE(hash(ab_1), hash(a));
  EXPECT_NE(hash(ab_1), hash(Vertex()));
}

//Tests for hash<std::pair>
TEST(HashPair, HashPair) {
  std::hash<std::pair<int, int> > hash;

  auto p11 = std::make_pair(1, 1);
  auto p00 = std::make_pair(0, 0);
  auto p01 = std::make_pair(0, 1);
  auto p10 = std::make_pair(1, 0);

  EXPECT_EQ(hash(p00), hash(std::make_pair(0, 0)));
  EXPECT_NE(hash(p00), hash(p11));
  EXPECT_NE(hash(p00), hash(p01));
  EXPECT_NE(hash(p00), hash(p10));
  EXPECT_NE(hash(p01), hash(p11));
  EXPECT_NE(hash(p01), hash(p10));
  EXPECT_NE(hash(p10), hash(p11));
}

}
}
}

