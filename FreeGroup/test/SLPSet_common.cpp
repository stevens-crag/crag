/*
 * SLPSet_common.cpp
 *
 *  Created on: Jan 21, 2013
 *      Author: dpantele
 */

#include "gtest/gtest.h"
#include "SLPSet.h"

#include <memory>
#include <functional>
#include <utility>

namespace crag{
namespace {
  class SignedVertexTest : public ::testing::Test {
    protected:
      SignedVertexTest()
        : v1p(SLPVertex::terminal_vertex(1))
        , v2p(SLPVertex::terminal_vertex(1))
        , null_vertex(SLPVertex::Null)
        , negative_null_vertex(SLPVertex::Null.negate())
      {
        v1n = v1p.negate();
      }

      SLPVertex v1p;
      SLPVertex v1n;
      SLPVertex v2p;
      SLPVertex null_vertex;
      SLPVertex negative_null_vertex;
  };

  TEST_F(SignedVertexTest, Equality) {
    EXPECT_EQ(v1p, v1p);
    EXPECT_NE(v1p, v2p);
    EXPECT_NE(v1p, v1n);
  }

  TEST_F(SignedVertexTest, NullVertex) {
    EXPECT_EQ(null_vertex, negative_null_vertex);

    EXPECT_NE(v1p, null_vertex);
    EXPECT_EQ(v1p.left_child(), null_vertex);
    EXPECT_EQ(v1p.left_child(), v1p.right_child());
    EXPECT_EQ(v1p.left_child(), v1p.right_child().negate());

    EXPECT_EQ(null_vertex.length(), 0);
    EXPECT_EQ(null_vertex.height(), 0);
    EXPECT_EQ(null_vertex.has_left_child(), false);
    EXPECT_EQ(null_vertex.has_right_child(), false);
    EXPECT_EQ(null_vertex.left_child(), null_vertex);
    EXPECT_EQ(null_vertex.right_child(), null_vertex);
    EXPECT_EQ(null_vertex.is_terminal(), false);
    EXPECT_EQ(null_vertex.terminal_symbol(), INVALID_TERMINAL);
  }

  TEST_F(SignedVertexTest, hash_function) {
    std::hash<SLPVertex> hash;
    EXPECT_EQ(hash(v1p), hash(v1p));
    EXPECT_NE(hash(v1p), hash(v1n));
    EXPECT_NE(hash(v1p), hash(v2p));

    EXPECT_EQ(hash(null_vertex), hash(negative_null_vertex));
    EXPECT_NE(hash(null_vertex), hash(v1p));
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

  TEST(SLPVertexTest, TerminalVertex) {
    SLPVertex a = SLPVertex::terminal_vertex(1);

    EXPECT_EQ(a.length(), 1);
    EXPECT_EQ(a.height(), 1);
    EXPECT_FALSE(a.has_left_child());
    EXPECT_FALSE(a.has_right_child());

    EXPECT_TRUE(a.is_terminal());
    EXPECT_EQ(a.terminal_symbol(), 1);
  }

  TEST(SLPVertexTest, ConcatednatedVertex) {
    SLPVertex a = SLPVertex::terminal_vertex(1);
    SLPVertex b = SLPVertex::terminal_vertex(2);
    SLPVertex c = SLPVertex::terminal_vertex(3);

    SLPVertex ab_1 = SLPVertex::concatenate(a, b.negate());
    EXPECT_EQ(ab_1.length(), 2);
    EXPECT_EQ(ab_1.height(), 2);
    EXPECT_EQ(ab_1.left_child(), a);
    EXPECT_NE(ab_1.right_child(), b);
    EXPECT_EQ(ab_1.right_child(), b.negate());

    SLPVertex d = SLPVertex::concatenate(ab_1.negate(), c.negate()); //ba^-1c_1
    EXPECT_EQ(d.length(), 3);
    EXPECT_EQ(d.height(), 3);
    EXPECT_EQ(d.left_child().left_child(), b);
    EXPECT_EQ(d.left_child().right_child(), a.negate());
    EXPECT_EQ(d.right_child(), c.negate());
  }

  TEST(LongInteger, Swap) {
    LongInteger a = 1;
    LongInteger b = 2;
    swap(a, b);
    EXPECT_EQ(a.get_si(), 2);
    EXPECT_EQ(b.get_si(), 1);
  }
}
}

