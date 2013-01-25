/*
 * SLPSet_inspector.cpp
 *
 *  Created on: Jan 22, 2013
 *      Author: dpantele
 */


#include "gtest/gtest.h"
#include "SLPSet.h"

#include <memory>
#include <functional>
#include <utility>

namespace crag {
namespace {
  TEST(Inspector, EmptyInspector) {
    EXPECT_TRUE(SLPPostorderInspector().inspection_ended());
  }

  TEST(Inspector, LeftTree) {
    SLPVertex a = SLPVertex::terminal_vertex(1);
    SLPVertex b = SLPVertex::terminal_vertex(2);
    SLPVertex c = SLPVertex::terminal_vertex(3);

    SLPVertex ab = SLPVertex::concatenate(a, b);
    SLPVertex abc = SLPVertex::concatenate(ab, c);

    SLPPostorderInspector inspector(abc);
    ASSERT_FALSE(inspector.inspection_ended());
    EXPECT_EQ(inspector.current_vertex(), a);
    inspector.go_to_next_vertex();
    ASSERT_FALSE(inspector.inspection_ended());
    EXPECT_EQ(inspector.current_vertex(), b);
    inspector.go_to_next_vertex();
    ASSERT_FALSE(inspector.inspection_ended());
    EXPECT_EQ(inspector.current_vertex(), ab);
    inspector.go_to_next_vertex();
    ASSERT_FALSE(inspector.inspection_ended());
    EXPECT_EQ(inspector.current_vertex(), c);
    inspector.go_to_next_vertex();
    ASSERT_FALSE(inspector.inspection_ended());
    EXPECT_EQ(inspector.current_vertex(), abc);
    inspector.go_to_next_vertex();
    EXPECT_TRUE(inspector.inspection_ended());
  }

  TEST(Inspector, RightTree) {
    SLPVertex a = SLPVertex::terminal_vertex(1);
    SLPVertex b = SLPVertex::terminal_vertex(2);
    SLPVertex c = SLPVertex::terminal_vertex(3);

    SLPVertex bc = SLPVertex::concatenate(b, c);
    SLPVertex abc = SLPVertex::concatenate(a, bc);

    SLPPostorderInspector inspector(abc);
    ASSERT_FALSE(inspector.inspection_ended());
    EXPECT_EQ(inspector.current_vertex(), a);
    inspector.go_to_next_vertex();
    ASSERT_FALSE(inspector.inspection_ended());
    EXPECT_EQ(inspector.current_vertex(), b);
    inspector.go_to_next_vertex();
    ASSERT_FALSE(inspector.inspection_ended());
    EXPECT_EQ(inspector.current_vertex(), c);
    inspector.go_to_next_vertex();
    ASSERT_FALSE(inspector.inspection_ended());
    EXPECT_EQ(inspector.current_vertex(), bc);
    inspector.go_to_next_vertex();
    ASSERT_FALSE(inspector.inspection_ended());
    EXPECT_EQ(inspector.current_vertex(), abc);
    inspector.go_to_next_vertex();
    EXPECT_TRUE(inspector.inspection_ended());
  }
  TEST(Inspector, CrossedTree) {
    SLPVertex a = SLPVertex::terminal_vertex(1);
    SLPVertex b = SLPVertex::terminal_vertex(2);

    SLPVertex ab = SLPVertex::concatenate(a, b);
    SLPVertex ba_1 = SLPVertex::concatenate(b, a.negate());//ba^-1

    SLPVertex abab_1 = SLPVertex::concatenate(ab, ba_1.negate());

    SLPPostorderInspector inspector(abab_1);
    ASSERT_FALSE(inspector.inspection_ended());
    EXPECT_EQ(inspector.current_vertex(), a);
    inspector.go_to_next_vertex();
    ASSERT_FALSE(inspector.inspection_ended());
    EXPECT_EQ(inspector.current_vertex(), b);
    inspector.go_to_next_vertex();
    ASSERT_FALSE(inspector.inspection_ended());
    EXPECT_EQ(inspector.current_vertex(), ab);
    inspector.go_to_next_vertex();
    ASSERT_FALSE(inspector.inspection_ended());
    EXPECT_EQ(inspector.current_vertex(), a);
    inspector.go_to_next_vertex();
    ASSERT_FALSE(inspector.inspection_ended());
    EXPECT_EQ(inspector.current_vertex(), b.negate());
    inspector.go_to_next_vertex();
    ASSERT_FALSE(inspector.inspection_ended());
    EXPECT_EQ(inspector.current_vertex(), ba_1.negate());
    inspector.go_to_next_vertex();
    ASSERT_FALSE(inspector.inspection_ended());
    EXPECT_EQ(inspector.current_vertex(), abab_1);
    inspector.go_to_next_vertex();
    EXPECT_TRUE(inspector.inspection_ended());
  }
}
}



