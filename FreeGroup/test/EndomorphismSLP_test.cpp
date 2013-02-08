/*
 * FreeGroupAutomorhpism_test.cpp
 *
 *  Created on: Jan 29, 2013
 *      Author: pmorar
 */

#include "gtest/gtest.h"
#include "../include/EndomorphismSLP.h"

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
template<int rank>
bool compareFGAutomorphisms(const EndomorphismSLP<rank>& a1, const EndomorphismSLP<rank>& a2) {
	if (&a1 == &a2)
			return true;
	for (TerminalSymbol i = 1; i < a1.generators_num(); ++i)
		if (!compare(a1.slp(i), a2.slp(i)))
			return false;
	return true;
}

class FreeGroupAutomorphismConstructorsTest : public ::testing::Test {
	static const int gen_num = 10;
    protected:
	FreeGroupAutomorphismConstructorsTest() :
		id(EndomorphismSLP<gen_num>::identity())
      {
      }
	EndomorphismSLP<gen_num> id;
//
//	bool is_nielsen1(const EndomorphismSLP& a) {
//		if (a.generators_num() != gen_num)
//					return false;
//		for (TerminalSymbol t = 1; t <= gen_num; ++t) {
//			if (t != changed_symbol && a.slp(t).terminal_symbol() != t)
//					return false;
//		}
//		SLPVertex v =  a.slp(changed_symbol);
//		return v.is_negative() && v.terminal_symbol() == changed_symbol;
//	}
//
//	bool is_nielsen2(const EndomorphismSLP& a) {
//		if (a.generators_num() != gen_num)
//			return false;
//		for (TerminalSymbol t = 1; t <= gen_num; ++t) {
//			if (t != changed_symbol && a.slp(t).terminal_symbol() != t)
//					return false;
//		}
//		SLPVertex v = a.slp(changed_symbol);
//		return v.has_left_child() && v.has_right_child() &&
//				v.left_child().is_terminal() &&
//				v.left_child().terminal_symbol() == left_child.terminal_symbol() &&
//				v.left_child().is_negative() == left_child.is_negative() &&
//				v.right_child().is_terminal() &&
//				v.right_child().terminal_symbol() == right_child.terminal_symbol() &&
//				v.right_child().is_negative() == right_child.is_negative();
//	}
//
  };
//
//TEST_F(FreeGroupAutomorphismConstructorsTest, SimpleConstructor) {
//for (int i = 0, n = 0; i < 10; ++i, n += i) {
//  EndomorphismSLP a(n);
//  ASSERT_EQ(a.generators_num(), n);
//}
//}
//
//TEST_F(FreeGroupAutomorphismConstructorsTest, NielsenConstructors) {
//	EXPECT_TRUE(is_nielsen1(nielsen1));
//	EXPECT_TRUE(is_nielsen2(nielsen2));
//}
//
//TEST_F(FreeGroupAutomorphismConstructorsTest, SimpleCopyConstructor) {
//	EndomorphismSLP a1(nielsen2);
//	EXPECT_EQ(a1.generators_num(), gen_num);
//	EXPECT_TRUE(is_nielsen2(a1));
//}
//
//TEST_F(FreeGroupAutomorphismConstructorsTest, SimpleMoveConstructor) {
//	EndomorphismSLP a1(std::move(nielsen2));
//	EXPECT_EQ(a1.generators_num(), gen_num);
//	EXPECT_EQ(nielsen2.generators_num(), 0);
//	EXPECT_TRUE(is_nielsen2(a1));
//}
//
//TEST_F(FreeGroupAutomorphismConstructorsTest, CompositionWithNielsen) {
//	EndomorphismSLP a1(gen_num);
//	a1.composeWithNielsen(changed_symbol);
//	EXPECT_TRUE(is_nielsen1(a1));
//
//	EndomorphismSLP a2(gen_num);
//	a2.composeWithNielsen(changed_symbol, left_child, right_child);
//	EXPECT_TRUE(is_nielsen2(a2));
//}
//
//
//
//TEST_F(FreeGroupAutomorphismConstructorsTest, ApplicationOfTwoNielsen) {
//	EndomorphismSLP a_ref(gen_num);
//	a_ref.composeWithNielsen(changed_symbol);
//	a_ref.composeWithNielsen(changed_symbol, left_child, right_child);
//
//	EndomorphismSLP a_ref2(a_ref);
//	a_ref2.composeWithNielsen(1, left_child, right_child);
//	a_ref2.composeWithNielsen(2, right_child, left_child);
//	EndomorphismSLP a_cp(a_ref2);
//	a_ref.composeWithNielsen(changed_symbol);
//	a_ref.composeWithNielsen(changed_symbol, left_child, right_child);
//
//	EndomorphismSLP a1(gen_num, changed_symbol, left_child, right_child);
//	a1.apply(EndomorphismSLP(gen_num, changed_symbol));
//	EXPECT_TRUE(compareFGAutomorphisms(a1, a_ref));
//
//	EndomorphismSLP a2(a_ref);
//	a2.apply(a_cp);
//	EXPECT_TRUE(compareFGAutomorphisms(a2, a_ref2));
//}
//




} /* namespace crag */
