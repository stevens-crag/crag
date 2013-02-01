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

} //namespace
} //namespace crag
