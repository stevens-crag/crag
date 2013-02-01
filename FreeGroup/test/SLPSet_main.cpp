/*
 * SLPSet_main.cpp
 *
 *  Created on: Jan 21, 2013
 *      Author: pmorar
 */

#include "gtest/gtest.h"
#include "SLPSet.h"


namespace crag {
class SLPSetConstructorsTest : public ::testing::Test {
    protected:
	SLPSetConstructorsTest()
        : v(SLPVertex::concatenate(
        		SLPVertex::terminal_vertex(1),
        		SLPVertex::concatenate(
        				SLPVertex::terminal_vertex(2),
        				SLPVertex::terminal_vertex(1)))),
       slp(v)
      {
      }

    SLPVertex v;
  	SLPSet slp;


  };

TEST_F(SLPSetConstructorsTest, SimpleConstructor) {
for (int i = 0, n = 0; i < 10; ++i, n += i) {
  SLPSet slp(n);
  ASSERT_EQ(slp.roots_num(), n);
  for (int j = 0; j < n; ++j) {
	SLPProducedWord word = slp.produced_word(j);
	EXPECT_EQ(word.size(), 1);
	const SLPVertex& v = word[0];
	EXPECT_TRUE(v.is_terminal());
	EXPECT_EQ(v.terminal_symbol(), j + 1);
  }
}
}

TEST_F(SLPSetConstructorsTest, SingleRootConstructor) {
	EXPECT_EQ(slp.roots_num(), 1);
	EXPECT_EQ(slp.root(0), v);
}

TEST_F(SLPSetConstructorsTest, SimpleCopyConstructor) {
	SLPSet slp1(slp);
	EXPECT_EQ(slp1.roots_num(), 1);
	EXPECT_EQ(slp.root(0), v);
}

TEST_F(SLPSetConstructorsTest, SimpleMoveConstructor) {
	SLPSet slp1(std::move(slp));
	EXPECT_EQ(slp1.roots_num(), 1);
	EXPECT_EQ(slp.roots_num(), 0);
	EXPECT_EQ(slp1.root(0), v);
}

TEST_F(SLPSetConstructorsTest, RangeConstructor) {
	std::vector<SLPVertex> range(10, v);
	SLPSet slp1(range.begin(), range.end());
	EXPECT_EQ(slp1.roots_num(), 10);
	for (int i = 0; i < 10; ++i)
		EXPECT_EQ(slp1.root(i), v) << "i = " << i;
}

TEST_F(SLPSetConstructorsTest, initializerListConstructor) {
	SLPSet slp1 = {v, v, v, v, v};
	EXPECT_EQ(slp1.roots_num(), 5);
	for (int i = 0; i < 5; ++i)
		EXPECT_EQ(slp1.root(i), v);
}

 
}
