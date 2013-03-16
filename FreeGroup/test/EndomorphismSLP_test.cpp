/*
 * FreeGroupAutomorhpism_test.cpp
 *
 *  Created on: Jan 29, 2013
 *      Author: pmorar
 */

#include "gtest/gtest.h"
#include "../include/EndomorphismSLP.h"

namespace crag {

typedef typename slp::TerminalVertexTemplate<int> TVertex;

//Auxiliary functions

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


  //! Finds images of terminal symbols under the composition of given morphisms. Each of the morphisms is either inverter, or left or right multiplier
  std::vector<std::string> find_image(const std::vector<EndomorphismSLP<int>>& morphisms) {
    unsigned int rank = 0;
    for (auto morph: morphisms) {
      const EndomorphismSLP<int>::size_type r = morph.max_non_trivial_image_symbol();
      if (r > rank)
        rank = r;
    }
    using std::vector;
    if (rank == 0)
      return vector<std::string>();
    vector<vector<int>> images;
    images.reserve(rank);
    //filling with id images
    for (int i = 1; i <= rank; ++i)
      images.push_back({i});
    //applying morphisms one by one
    for (auto morph: morphisms) {
      vector<vector<int>> tmp;//we fill it with images under morph
      tmp.reserve(rank);
      for (int terminal = 1; terminal <= images.size(); ++terminal) {
        const slp::Vertex v = morph.slp(terminal);
        assert (v.height() <= 2);

        if (v.height() == 1) {//inverter or id
          const vector<int>& img = images[terminal - 1];
          const int terminal_symbol = TVertex(v).terminal_symbol();
          if (terminal_symbol == terminal) {//id
            tmp.push_back(img);
          } else {//inverter
            assert(-terminal_symbol == terminal);
            vector<int> tmp_img;
            tmp_img.reserve(img.size());
            for (auto r_it = img.rbegin(); r_it != img.rend(); ++r_it) {
              tmp_img.push_back(-*r_it);
            }
            tmp.push_back(tmp_img);
          }
        } else {//righ_or left multiplier
          const int left_terminal_symbol = TVertex(v.left_child()).terminal_symbol();
          const int right_terminal_symbol = TVertex(v.right_child()).terminal_symbol();
          assert(left_terminal_symbol > 0 && right_terminal_symbol > 0);

          vector<int> l_img = left_terminal_symbol <= rank ? images[left_terminal_symbol - 1] : vector<int>({left_terminal_symbol});
          vector<int> r_img = right_terminal_symbol <= rank ? images[right_terminal_symbol - 1] : vector<int>({right_terminal_symbol});
          vector<int> tmp_img;
          tmp_img.reserve(l_img.size() + r_img.size());
          tmp_img.insert(tmp_img.end(), l_img.begin(), l_img.end());
          tmp_img.insert(tmp_img.end(), r_img.begin(), r_img.end());
          tmp.push_back(tmp_img);
        }
      }
      using std::swap;
      swap(tmp, images);
    }

    std::vector<std::string> s_images;
    s_images.reserve(rank);

    for (auto v_img: images) {
      std::stringstream s_str;
      for (auto i: v_img)
        s_str << i;
      s_images.push_back(s_str.str());
    }

    return s_images;
  }

class EndomorphismSLPTest : public ::testing::Test {
    protected:
  typedef EndomorphismSLP<int> EMorphism;

  class CachedProducer {
  public:
   CachedProducer(const std::vector<EMorphism>& morph)
    : morphisms_(morph), i(0) {}

   EMorphism operator()() {
     return morphisms_[i++];
   }

  private:
   const std::vector<EMorphism>& morphisms_;
   int i;
  };
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

TEST_F(EndomorphismSLPTest, DirectImageCalculationTest) {
  EXPECT_EQ(0, find_image({}).size());
  EXPECT_EQ(0, find_image({EMorphism::identity()}).size());

  std::vector<EMorphism> morphisms = {EMorphism::right_multiplier(1, 2)};
  auto img = find_image(morphisms);
  EXPECT_EQ(1, img.size());
  EXPECT_EQ("12", img[0]);

  morphisms.push_back(EMorphism::right_multiplier(1, 2));
  img = find_image(morphisms);
  EXPECT_EQ(1, img.size());
  EXPECT_EQ("122", img[0]);

  morphisms.push_back(EMorphism::inverter(3));
  img = find_image(morphisms);
  EXPECT_EQ(3, img.size());
  EXPECT_EQ("122", img[0]);
  EXPECT_EQ("2", img[1]);
  EXPECT_EQ("-3", img[2]);

  morphisms.push_back(EMorphism::left_multiplier(3, 1));
  img = find_image(morphisms);
  EXPECT_EQ(3, img.size());
  EXPECT_EQ("-3122", img[0]);
  EXPECT_EQ("2", img[1]);
  EXPECT_EQ("-3", img[2]);

  morphisms.push_back(EMorphism::left_multiplier(1, 2));
  img = find_image(morphisms);
  EXPECT_EQ(3, img.size());
  EXPECT_EQ("-3122", img[0]);
  EXPECT_EQ("-31222", img[1]);
  EXPECT_EQ("-3", img[2]);

  morphisms.push_back(EMorphism::inverter(2));
  img = find_image(morphisms);
  EXPECT_EQ(3, img.size());
  EXPECT_EQ("-3122", img[0]);
  EXPECT_EQ("-2-2-2-13", img[1]);
  EXPECT_EQ("-3", img[2]);
}

TEST_F(EndomorphismSLPTest, RandomGeneratorConstructorsTest) {
  UniformAutomorphismSLPGenerator<int> rnd(5);//default
  UniformAutomorphismSLPGenerator<int> rnd1(10, 546457);//with seed
  std::default_random_engine r_engine;
  UniformAutomorphismSLPGenerator<int> rnd2(4, &r_engine);//with generator
  for (int i = 0; i < 100; ++i) {
      EXPECT_EQ(1, rnd().non_trivial_images_num());
      EXPECT_EQ(1, rnd1().non_trivial_images_num());
      EXPECT_EQ(1, rnd2().non_trivial_images_num());
    }
}

TEST_F(EndomorphismSLPTest, RandomGeneratorStressTest) {
  for (auto rank : {1, 5, 10, 20, 30}) {
    UniformAutomorphismSLPGenerator<int> rnd(rank);
    for (auto size : {0, 5, 10, 30, 60}) {
      std::vector<EMorphism> morphisms;
      for (int i = 0; i < size; ++i)
        morphisms.push_back(rnd());

      auto morphism = EMorphism::composition(morphisms.begin(), morphisms.end());
      auto direct_images = find_image(morphisms);
      ASSERT_GE(rank, direct_images.size());
      ASSERT_GE(rank, morphism.max_non_trivial_image_symbol());
      for (int i = 0; i < direct_images.size(); ++i)
        EXPECT_EQ(direct_images[i], to_string(morphism.image(i + 1)));
    }
  }
}


TEST_F(EndomorphismSLPTest, ProducerVsIteratorGenerationEquality) {
  for (auto rank : {1, 5, 10, 20, 30}) {
      UniformAutomorphismSLPGenerator<int> rnd(rank);
      for (auto size : {0, 5, 10, 30, 60}) {
      std::vector<EMorphism> morphisms;
      for (int i = 0; i < size; ++i)
        morphisms.push_back(rnd());

      auto morphism1 = EMorphism::composition(morphisms.begin(), morphisms.end());
      CachedProducer cp(morphisms);
      auto morphism2 = EMorphism::composition(size, cp);
      for (int i = 0; i < rank; ++i)
        EXPECT_EQ(to_string(morphism1.image(i + 1)), to_string(morphism2.image(i + 1)));
    }
  }
}


} /* namespace crag */
