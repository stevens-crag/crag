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

namespace {
  TEST(Inspector, EmptyInspector) {
    EXPECT_TRUE(SLPPostorderInspector().inspection_ended());
  }

  TEST(Inspector, SomeTreeInspection) {

  }
}




