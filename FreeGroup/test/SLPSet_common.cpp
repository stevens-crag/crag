/*
 * SLPSet_common.cpp
 *
 *  Created on: Jan 21, 2013
 *      Author: dpantele
 */

#include "gtest/gtest.h"
#include "SLPSet.cpp"

#include <memory>
#include <functional>
#include <utility>

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
}


