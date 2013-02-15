/*
 * SLPSet_common.cpp
 *
 *  Created on: Jan 21, 2013
 *      Author: dpantele
 */

#include "gtest/gtest.h"
#include "slp.h"

#include <memory>
#include <functional>
#include <utility>

namespace crag {
namespace slp {
namespace {
  typedef TerminalVertexTemplate<int> TerminalVertex;
  class SlpVertexTest : public ::testing::Test {
    protected:
      SlpVertexTest()
        : v1p(TerminalVertex::create(1))
        , v1n(TerminalVertex::create(-1))
        , v2p(TerminalVertex::create(2))
        , null_vertex(Vertex::Null)
        , negative_null_vertex(Vertex::Null->negate())
      { }

      Vertex::Ptr v1p;
      Vertex::Ptr v1n;
      Vertex::Ptr v2p;
      Vertex::Ptr null_vertex;
      Vertex::Ptr negative_null_vertex;
  };

  TEST_F(SlpVertexTest, Equality) {
    EXPECT_EQ(*v1p, *v1p);
    EXPECT_NE(*v1p, *v2p);
    EXPECT_NE(*v1p, *v1n);
  }

  TEST_F(SlpVertexTest, NullVertex) {
    EXPECT_EQ(*null_vertex, *negative_null_vertex);

    EXPECT_NE(*v1p, *null_vertex);
    EXPECT_EQ(*v1p->left_child(), *null_vertex);
    EXPECT_EQ(*v1p->left_child(), *v1p->right_child());
    EXPECT_EQ(*v1p->left_child(), *v1p->right_child()->negate());

    EXPECT_EQ(null_vertex->length(), 0);
    EXPECT_EQ(null_vertex->height(), 0);
    EXPECT_EQ(null_vertex->has_left_child(), false);
    EXPECT_EQ(null_vertex->has_right_child(), false);
    EXPECT_EQ(*null_vertex->left_child(), *null_vertex);
    EXPECT_EQ(*null_vertex->right_child(), *null_vertex);
  }

  TEST_F(SlpVertexTest, hash_function) {
    ::std::hash<Vertex> hash;
    EXPECT_EQ(hash(*v1p), hash(*v1p));
    EXPECT_NE(hash(*v1p), hash(*v1n));
    EXPECT_NE(hash(*v1p), hash(*v2p));

    EXPECT_EQ(hash(*null_vertex), hash(*negative_null_vertex));
    EXPECT_NE(hash(*null_vertex), hash(*v1p));
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

  TEST(SlpVertexTypesTest, TerminalVertexTest) {
    Vertex::Ptr a = TerminalVertex::create(1);

    EXPECT_EQ(a->length(), 1);
    EXPECT_EQ(a->height(), 1);
    EXPECT_FALSE(a->has_left_child());
    EXPECT_FALSE(a->has_right_child());

    //EXPECT_EQ(1, *std::dynamic_pointer_cast<TerminalVertex>(a));
  }

  TEST(SlpVertexTypesTest, ConcatednatedVertex) {
    auto a = TerminalVertex::create(1);
    auto b = TerminalVertex::create(2);
    auto c = TerminalVertex::create(3);

    auto ab_1 = NonterminalVertex::create(a, b->negate());
    EXPECT_EQ(ab_1->length(), 2);
    EXPECT_EQ(ab_1->height(), 2);
    EXPECT_EQ(*ab_1->left_child(), *a);
    EXPECT_NE(*ab_1->right_child(), *b);
    EXPECT_EQ(*ab_1->right_child(), *b->negate());

    auto d = NonterminalVertex::create(ab_1->negate(), c->negate()); //ba^-1c_1
    EXPECT_EQ(d->length(), 3);
    EXPECT_EQ(d->height(), 3);
    EXPECT_EQ(*d->left_child()->negate(), *ab_1);
    EXPECT_EQ(*d->left_child()->left_child(), *b);
    EXPECT_EQ(*d->left_child()->right_child(), *a->negate());
    EXPECT_EQ(*d->right_child(), *c->negate());
  }
}
}
}

