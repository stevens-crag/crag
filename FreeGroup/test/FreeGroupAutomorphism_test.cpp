/*
 * FreeGroupAutomorhpism_test.cpp
 *
 *  Created on: Jan 29, 2013
 *      Author: pmorar
 */

#include "gtest/gtest.h"
#include "../include/FreeGroupAutomorhpism.h"

namespace crag {


class FreeGroupAutomorphismConstructorsTest : public ::testing::Test {
    protected:
	FreeGroupAutomorphismConstructorsTest() :
		gen_num(5),
		changed_symbol(3),
		left_child(SLPVertex::terminal_vertex(2).negate()),
		right_child(SLPVertex::terminal_vertex(1)),
		nielsen1(gen_num, changed_symbol),
		nielsen2(gen_num, changed_symbol, left_child, right_child)
      {
      }

	FreeGroupAutomorphism::size_type gen_num;
	TerminalSymbol changed_symbol;
	SLPVertex left_child, right_child;
	FreeGroupAutomorphism nielsen1;
	FreeGroupAutomorphism nielsen2;

	bool is_nielsen1(const FreeGroupAutomorphism& a) {
		if (a.generators_num() != gen_num)
					return false;
		for (TerminalSymbol t = 1; t <= gen_num; ++t) {
			if (t != changed_symbol && a.slp(t).terminal_symbol() != t)
					return false;
		}
		SLPVertex v =  a.slp(changed_symbol);
		return v.is_negative() && v.terminal_symbol() == changed_symbol;
	}

	bool is_nielsen2(const FreeGroupAutomorphism& a) {
		if (a.generators_num() != gen_num)
			return false;
		for (TerminalSymbol t = 1; t <= gen_num; ++t) {
			if (t != changed_symbol && a.slp(t).terminal_symbol() != t)
					return false;
		}
		SLPVertex v =  a.slp(changed_symbol);
		return v.has_left_child() && v.has_right_child() &&
				v.left_child() == left_child && v.right_child() == right_child;
	}

  };

TEST_F(FreeGroupAutomorphismConstructorsTest, SimpleConstructor) {
for (int i = 0, n = 0; i < 10; ++i, n += i) {
  FreeGroupAutomorphism a(n);
  ASSERT_EQ(a.generators_num(), n);
}
}

TEST_F(FreeGroupAutomorphismConstructorsTest, NielsenConstructors) {
	EXPECT_TRUE(is_nielsen1(nielsen1));
	EXPECT_TRUE(is_nielsen2(nielsen2));
}

TEST_F(FreeGroupAutomorphismConstructorsTest, SimpleCopyConstructor) {
	FreeGroupAutomorphism a1(nielsen2);
	EXPECT_EQ(a1.generators_num(), gen_num);
	EXPECT_TRUE(is_nielsen2(a1));
}

TEST_F(FreeGroupAutomorphismConstructorsTest, SimpleMoveConstructor) {
	FreeGroupAutomorphism a1(std::move(nielsen2));
	EXPECT_EQ(a1.generators_num(), gen_num);
	EXPECT_EQ(nielsen2.generators_num(), 0);
	EXPECT_TRUE(is_nielsen2(a1));
}




} /* namespace crag */
