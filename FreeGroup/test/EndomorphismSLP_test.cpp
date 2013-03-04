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
bool compare_endomorphisms_directly(const EndomorphismSLP<int>& e1, const EndomorphismSLP<int>& e2) {
	if (&e1 == &e2)
			return true;
	const int max = e1.max_non_trivial_image_symbol();
	if (max != e2.max_non_trivial_image_symbol())
	  return false;
	for (int ts = 1; ts < max; ++ts) {
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
std::ostream& operator<<(std::ostream& out, const slp::VertexWord<int>& word) {
  for (auto wi = word.begin(); wi != word.end(); ++wi) {
    out << *wi;
  }
  return out;
}

//! Return string representation of vertex
std::string to_string(const slp::VertexWord<int>& word) {
  std::stringstream out;
  out << word;
  return out.str();
}

//! Prints images of terminal symbols
std::ostream& operator<<(std::ostream& out, const EndomorphismSLP<int>& e) {
  const int max = e.max_non_trivial_image_symbol();
  out << "[max non trivial symbol=" << max;
  for (int ts = 1; ts <= max; ++ts)
    out << "," << std::endl << ts << " => " << e.image(ts);
  return out << "]" << std::endl;
}

////! Prints images of terminal symbols
//std::string find_image(unsigned int symbol, const std::vector<EndomorphismSLP<int>>& morphisms) {
//  std::string intermediate_image = symbol;
//  const int max = e.max_non_trivial_image_symbol();
//  out << "[max non trivial symbol=" << max;
//  for (int ts = 1; ts <= max; ++ts)
//    out << "," << std::endl << ts << " => " << e.image(ts);
//  return out << "]" << std::endl;
//}


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
  auto prod = EMorphism::right_multiplier(1, 2) * EMorphism::identity();
  EXPECT_EQ(1, prod.non_trivial_images_num()) << prod;
  EXPECT_EQ("12", to_string(prod.image(1))) << prod;

  prod = EMorphism::identity() * EMorphism::left_multiplier(1, 2);
  EXPECT_EQ(1, prod.non_trivial_images_num()) << prod;
  EXPECT_EQ("12", to_string(prod.image(2))) << prod;

  prod = EMorphism::right_multiplier(1, 2) * EMorphism::right_multiplier(1, 2);
  EXPECT_EQ(1, prod.non_trivial_images_num()) << prod;
  EXPECT_EQ("122", to_string(prod.image(1))) << prod;

  prod = EMorphism::right_multiplier(1, 2) * EMorphism::left_multiplier(1, 2);
  EXPECT_EQ(2, prod.non_trivial_images_num()) << prod;
  EXPECT_EQ("12", to_string(prod.image(1))) << prod;
  EXPECT_EQ("122", to_string(prod.image(2))) << prod;

  prod = EMorphism::right_multiplier(1, 2) * EMorphism::inverter(1);
  EXPECT_EQ(1, prod.non_trivial_images_num()) << prod;
  EXPECT_EQ("-2-1", to_string(prod.image(1))) << prod;

  prod = EMorphism::inverter(1) * EMorphism::left_multiplier(1, 2);
  EXPECT_EQ(2, prod.non_trivial_images_num()) << prod;
  EXPECT_EQ("-1", to_string(prod.image(1))) << prod;
  EXPECT_EQ("-12", to_string(prod.image(2))) << prod;
}


TEST_F(EndomorphismSLPTest, Composition1) {
  auto e = EMorphism::right_multiplier(1, 2);
  e *= EMorphism::right_multiplier(1, 2);
  EXPECT_EQ(1, e.non_trivial_images_num());
  EXPECT_EQ("122", to_string(e.image(1))) << e;

  e *= EMorphism::inverter(3);
  EXPECT_EQ(2, e.non_trivial_images_num());
  EXPECT_EQ("122", to_string(e.image(1))) << e;
  EXPECT_EQ("-3", to_string(e.image(3))) << e;

  e *= EMorphism::left_multiplier(3, 1);
  EXPECT_EQ(2, e.non_trivial_images_num());
  EXPECT_EQ("-3122", to_string(e.image(1))) << e;
  EXPECT_EQ("-3", to_string(e.image(3))) << e;

  e *= EMorphism::left_multiplier(1, 2);
  EXPECT_EQ(3, e.non_trivial_images_num());
  EXPECT_EQ("-3122", to_string(e.image(1))) << e;
  EXPECT_EQ("-31222", to_string(e.image(2))) << e;
  EXPECT_EQ("-3", to_string(e.image(3))) << e;

  e *= EMorphism::inverter(2);
  EXPECT_EQ(3, e.non_trivial_images_num());
  EXPECT_EQ("-3122", to_string(e.image(1))) << e;
  EXPECT_EQ("-2-2-2-13", to_string(e.image(2))) << e;
  EXPECT_EQ("-3", to_string(e.image(3))) << e;
}


} /* namespace crag */
