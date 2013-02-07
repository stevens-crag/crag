/*
 * FreeGroupAutomorhpism_test.cpp
 *
 *  Created on: Jan 29, 2013
 *      Author: pmorar
 */

#include "gtest/gtest.h"
#include "../include/FreeGroupAutomorhpism.h"

namespace crag {

//! Compares two vertices by comparing two words produced by them
bool compare(const SLPVertex& v1, const SLPVertex& v2) {
	if (&v1 == &v2)
		return true;
	if (v1.length() != v2.length())
		return false;
	SLPProducedWord w1(v1);
	SLPProducedWord w2(v2);
	for (auto it1 = w1.begin(), it2 = w2.begin();
			it1 != w1.end();//already checked that the length are the same
			++it1, ++it2) {
		auto t1 = *it1;
		auto t2 = *it2;
		if (t1.is_negative() != t2.is_negative() || t1.terminal_symbol() != t2.terminal_symbol())
			return false;
	}
	return true;
}

//! Compares two automorphisms
bool compareFGAutomorphisms(const FreeGroupAutomorphism& a1, const FreeGroupAutomorphism& a2) {
	if (&a1 == &a2)
			return true;
	if (a1.generators_num() != a2.generators_num())
		return false;
	for (TerminalSymbol i = 1; i < a1.generators_num(); ++i)
		if (!compare(a1.slp(i), a2.slp(i)))
			return false;
	return true;
}

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
		SLPVertex v = a.slp(changed_symbol);
		return v.has_left_child() && v.has_right_child() &&
				v.left_child().is_terminal() &&
				v.left_child().terminal_symbol() == left_child.terminal_symbol() &&
				v.left_child().is_negative() == left_child.is_negative() &&
				v.right_child().is_terminal() &&
				v.right_child().terminal_symbol() == right_child.terminal_symbol() &&
				v.right_child().is_negative() == right_child.is_negative();
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

TEST_F(FreeGroupAutomorphismConstructorsTest, CompositionWithNielsen) {
	FreeGroupAutomorphism a1(gen_num);
	a1.composeWithNielsen(changed_symbol);
	EXPECT_TRUE(is_nielsen1(a1));

	FreeGroupAutomorphism a2(gen_num);
	a2.composeWithNielsen(changed_symbol, left_child, right_child);
	EXPECT_TRUE(is_nielsen2(a2));
}



TEST_F(FreeGroupAutomorphismConstructorsTest, ApplicationOfTwoNielsen) {
	FreeGroupAutomorphism a_ref(gen_num);
	a_ref.composeWithNielsen(changed_symbol);
	a_ref.composeWithNielsen(changed_symbol, left_child, right_child);

	FreeGroupAutomorphism a_ref2(a_ref);
	a_ref2.composeWithNielsen(1, left_child, right_child);
	a_ref2.composeWithNielsen(2, right_child, left_child);
	FreeGroupAutomorphism a_cp(a_ref2);
	a_ref.composeWithNielsen(changed_symbol);
	a_ref.composeWithNielsen(changed_symbol, left_child, right_child);

	FreeGroupAutomorphism a1(gen_num, changed_symbol, left_child, right_child);
	a1.apply(FreeGroupAutomorphism(gen_num, changed_symbol));
	EXPECT_TRUE(compareFGAutomorphisms(a1, a_ref));

	FreeGroupAutomorphism a2(a_ref);
	a2.apply(a_cp);
	EXPECT_TRUE(compareFGAutomorphisms(a2, a_ref2));
}





} /* namespace crag */
