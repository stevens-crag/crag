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
      virtual void SetUp() {
        auto vertex1 = std::make_shared<SLPVertex>();
        auto vertex2 = std::make_shared<SLPVertex>();

        v1p.ptr = vertex1;
        v1p.negative = false;

        v1n.ptr = vertex1;
        v1n.negative = true;

        v2p.ptr = vertex2;
        v2p.negative = false;

        null_vertex = SignedVertex::Null;
        negative_null_vertex = SignedVertex::Null;
        negative_null_vertex.negative = true;
      }

      SignedVertex v1p;
      SignedVertex v1n;
      SignedVertex v2p;
      SignedVertex null_vertex;
      SignedVertex negative_null_vertex;
  };

  TEST_F(SignedVertexTest, Equality) {
    EXPECT_EQ(v1p, v1p);
    EXPECT_NE(v1p, v2p);
    EXPECT_NE(v1p, v1n);
  }

  TEST_F(SignedVertexTest, NullVertex) {
    EXPECT_EQ(null_vertex, negative_null_vertex);

    EXPECT_NE(v1p, null_vertex);
    EXPECT_EQ(v1p.ptr->left_child(), null_vertex);
    EXPECT_EQ(v1p.ptr->left_child(), v1p.ptr->right_child());

  }

  TEST_F(SignedVertexTest, hash_function) {
    std::hash<SignedVertex> hash;
    EXPECT_EQ(hash(v1p), hash(v1p));
    EXPECT_NE(hash(v1p), hash(v1n));
    EXPECT_NE(hash(v1p), hash(v2p));

    EXPECT_EQ(hash(null_vertex), hash(negative_null_vertex));
    EXPECT_NE(hash(null_vertex), hash(v1p));
  }
}


