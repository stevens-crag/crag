/*
 * FreeGroupAutomorhpism_test.cpp
 *
 *  Created on: Jan 29, 2013
 *      Author: pmorar
 */

#include "gtest/gtest.h"
#include "../include/EndomorphismSLP.h"

namespace crag {


//! Compares two endomorphisms by comparing images as strings -- very inefficient
template<typename TerminalSymbol>
bool compare_endomorphisms_directly(const EndomorphismSLP<TerminalSymbol>& e1, const EndomorphismSLP<TerminalSymbol>& e2) {
	if (&e1 == &e2)
			return true;
	const TerminalSymbol max = e1.max_non_trivial_image_symbol();
	if (max != e2.max_non_trivial_image_symbol())
	  return false;
	for (TerminalSymbol ts = TerminalSymbol(1); ts < max; ++ts) {
	  auto word1 = e1.image(ts);
	  auto word2 = e2.image(ts);
	  if (word1.size() != word2.size())
	    return false;
	  for (auto wi1 = word1.begin(), wi2 = word2.begin();
	      wi1 != word1.end() && wi2 != word2.end();
	      ++wi1, ++wi2) {
	    if (*wi1 != *wi2)
	      return false;
	  }
	}
	return true;
}

//! Prints vertex
template<typename TerminalSymbol>
std::ostream& operator<<(std::ostream& out, const slp::VertexWord<TerminalSymbol>& word) {
  for (auto wi = word.begin(); wi != word.end(); ++wi) {
    out << *wi;
  }
  return out;
}

//! Return string representation of vertex
template<typename TerminalSymbol>
std::string to_string(const slp::VertexWord<TerminalSymbol>& word) {
  std::stringstream out;
  out << word;
  return out.str();
}

//! Prints images of terminal symbols
template<typename TerminalSymbol>
std::ostream& operator<<(std::ostream& out, const EndomorphismSLP<TerminalSymbol>& e) {
  const TerminalSymbol max = e.max_non_trivial_image_symbol();
  out << "[max non trivial symbol=" << max;
  for (TerminalSymbol ts = TerminalSymbol(1); ts < max; ++ts)
    out << "," << std::endl << e.image(ts);
  return out << "]";
}


class EndomorphismSLPTest : public ::testing::Test {
    protected:
  typedef EndomorphismSLP<int> EMorphism;
};

TEST_F(EndomorphismSLPTest, Identity) {
  auto id = EMorphism::identity();
  ASSERT_EQ(id.max_non_trivial_image_symbol(), int(0));
  EXPECT_EQ("1", to_string(id.image(1)));
  EXPECT_EQ("-1", to_string(id.image(-1)));
}

TEST_F(EndomorphismSLPTest, Inverter) {
  for (int i = 1; i < 10; ++i) {
    auto inverter = EMorphism::inverter(i);
    EXPECT_EQ(1, inverter.non_trivial_images_num());
    auto img = inverter.image(i);
    EXPECT_EQ(1, img.size());
    EXPECT_EQ(-i, img[0]);
  }
}

TEST_F(EndomorphismSLPTest, LeftMultiplier) {
  for (int i = 1; i < 10; ++i)
    for (int j = 1; j < 10; ++j) {
      if (i == j) 
        continue;
      auto inverter = EMorphism::left_multiplier(i, j);
      EXPECT_EQ(1, inverter.non_trivial_images_num());
      auto img = inverter.image(j);
      EXPECT_EQ(2, img.size());
      EXPECT_EQ(i, img[0]);
      EXPECT_EQ(j, img[1]);
  }
}

TEST_F(EndomorphismSLPTest, RightMultiplier) {
  for (int i = 1; i < 10; ++i)
    for (int j = 1; j < 10; ++j) {
      if (i == j)
        continue;
      auto inverter = EMorphism::right_multiplier(i, j);
      EXPECT_EQ(1, inverter.non_trivial_images_num());
      auto img = inverter.image(i);
      EXPECT_EQ(2, img.size());
      EXPECT_EQ(i, img[0]);
      EXPECT_EQ(j, img[1]);
  }
}

TEST_F(EndomorphismSLPTest, BasicComposition) {
  auto e = EMorphism::right_multiplier(1, 2);
  std::cout << "checking identites composition " << std::endl;
  std::cout << EMorphism::identity();
  std::cout << "checking other stuff" << std::endl;
  auto prod = EMorphism::right_multiplier(1, 2) * EMorphism::identity();
//  auto prod = EMorphism::identity() * EMorphism::right_multiplier(1, 2);

//  std::cout << to_string(prod.image(2)) << std::endl;
//  std::cout << to_string(prod.image(1)) << std::endl;
//  EXPECT_EQ(1, prod.max_non_trivial_image_symbol());
//  EXPECT_EQ("12", to_string(prod.image(1)));
//  EXPECT_TRUE(compare_endomorphisms_directly(e, prod)) << "e = " << e << ", prod = " << prod;
//  prod = EMorphism::right_multiplier(1, 2) * EMorphism::identity();
//  EXPECT_TRUE(compare_endomorphisms_directly(e, prod))  << "e = " << e << ", prod = " << prod;
}


//TEST_F(EndomorphismSLPTest, BasicCompositionCommutativity) {
//    std::vector<EMorphism> operand1, operand2;
//    //filling operands
//    for (int i = 1; i < 5; ++i)
//        for (int j = 1; j < 5; ++j) {
//          if (i == j)
//            continue;
//          auto right = EMorphism::right_multiplier(i, j);
//          auto left = EMorphism::left_multiplier(i, j);
//          operand1.push_back(left);
//          operand1.push_back(right);
//          operand2.push_back(left);
//          operand2.push_back(right);
//      }
//    //checking commutativity
//    for (auto e1: operand1)
//      for (auto e2: operand2) {
//        auto prod1 = e1 * e2;
//        auto prod2 = e2 * e1;
//        EXPECT_TRUE(compare_endomorphisms_directly(e1 * e2, e2 * e1))
//          << "e1*e2=" << prod1 << std::endl << "e2*e1=" << prod2 << std::endl;
//      }
//}


TEST_F(EndomorphismSLPTest, Composition1) {
  auto e = EMorphism::right_multiplier(1, 2);
  e *= EMorphism::right_multiplier(1, 2);
  EXPECT_EQ(1, e.non_trivial_images_num());
  EXPECT_EQ("122", to_string(e.image(1))) << "Automorphism = " << e;

//  e *= EMorphism::inverter(3);
//  EXPECT_EQ(2, e.non_trivial_images_num());
//  EXPECT_EQ("122", to_string(e.image(1))) << "Automorphism = " << e;
//  EXPECT_EQ("-3", to_string(e.image(3))) << "Automorphism = " << e;
//
//  e *= EMorphism::left_multiplier(3, 1);
//  EXPECT_EQ(2, e.non_trivial_images_num());
//  EXPECT_EQ("-3122", to_string(e.image(1))) << "Automorphism = " << e;
//  EXPECT_EQ("-3", to_string(e.image(3))) << "Automorphism = " << e;
//
//  e *= EMorphism::left_multiplier(1, 2);
//  EXPECT_EQ(3, e.non_trivial_images_num());
//  EXPECT_EQ("-3122", to_string(e.image(1))) << "Automorphism = " << e;
//  EXPECT_EQ("-3122", to_string(e.image(2))) << "Automorphism = " << e;
//  EXPECT_EQ("-3", to_string(e.image(3))) << "Automorphism = " << e;
}


} /* namespace crag */
