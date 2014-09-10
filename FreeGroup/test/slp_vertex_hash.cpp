/*
 * SLP_vertex_word.cpp
 *
 *  Created on: Feb 21, 2013
 *      Author: dpantele
 */

#include <memory>
#include <functional>
#include <utility>
#include <algorithm>
#include <string>

#include "Permutation.h"
#include "gtest/gtest.h"
#include "slp_vertex_hash.h"
#include "EndomorphismSLP.h"

namespace crag {
namespace slp {
namespace {

template <typename T>
class HasherTest : public ::testing::Test {
  protected:
    typedef slp::TerminalVertex TerminalVertex;
    TerminalVertex a;
    TerminalVertex b;

    HasherTest()
      : a(1)
      , b(2)
    { }

    virtual ~HasherTest() {}
};


TYPED_TEST_CASE_P(HasherTest);

TYPED_TEST_P(HasherTest, InterfaceImplementation) {
  typedef TypeParam VertexHash;

  typename HasherTest<TypeParam>::TerminalVertex& a = this->a;
  typename HasherTest<TypeParam>::TerminalVertex& b = this->b;

  EXPECT_TRUE(VertexHash(a.vertex_id()).is_equal_to(VertexHash(a.vertex_id())));
  EXPECT_TRUE(!VertexHash(a.vertex_id()).is_equal_to(VertexHash(b.vertex_id())));

  EXPECT_TRUE(VertexHash().is_equal_to(VertexHash()));

  VertexHash a_hash = VertexHash(a.vertex_id());
  EXPECT_TRUE(a_hash.is_equal_to(VertexHash(a_hash)));

  VertexHash ba = VertexHash(b.vertex_id());
  ba.concatenate_with(a_hash);
  EXPECT_TRUE(!VertexHash(a.vertex_id()).is_equal_to(ba));

  VertexHash a_inverse = VertexHash(a_hash);
  a_inverse.inverse_inplace();
  EXPECT_TRUE(!a_hash.is_equal_to(a_inverse));
  EXPECT_TRUE(a_inverse.is_equal_to(VertexHash(-a.vertex_id())));
}

REGISTER_TYPED_TEST_CASE_P(HasherTest, InterfaceImplementation);

typedef ::testing::Types<hashers::SinglePowerHash, hashers::PermutationHash<Permutation16> > HasherTypes;
INSTANTIATE_TYPED_TEST_CASE_P(All, HasherTest, HasherTypes);

TEST(HasherTest, SimpleCombinedHash) {
  typedef TVertexHash<hashers::SinglePowerHash, hashers::PermutationHash<Permutation16> > VertexHash;
  auto a = VertexHash(1);
  auto b = VertexHash(2);

  EXPECT_NE(a, b);

  auto ab = VertexHash(a) *= b;
  auto ba = VertexHash(b) *= a;

  EXPECT_NE(ab, ba);

  auto b_1a_1 = (VertexHash(-2) *= VertexHash(-1));

  EXPECT_EQ(b_1a_1, ab.inverse());
}

std::string print_tree_preorder(const Vertex& vertex) {
  std::ostringstream out;
  PreorderInspector inspector(vertex);
  while (!inspector.stopped()) {
    PrintTo(inspector.vertex(), &out);
    out << " -- ";

    ++inspector;
  }

  return out.str();
}

std::string print_tree_preorder_single(const Vertex& vertex) {
  std::ostringstream out;
  std::unordered_map<slp::Vertex, bool> mapping;
  auto acceptor = [&] (const inspector::InspectorTask& task) {
    return mapping.find(task.vertex) == mapping.end();//do not accept if vertex is mapped already
  };

  Inspector<inspector::Preorder, decltype(acceptor)> inspector(vertex, acceptor);
  while (!inspector.stopped()) {
    PrintTo(inspector.vertex(), &out);
    out << " << ";
    if (mapping.count(inspector.vertex().left_child()) || mapping.count(inspector.vertex().left_child().negate())) {
      out << "(p) ";
    }
    PrintTo(inspector.vertex().left_child(), &out);
    out << " >> ";
    if (mapping.count(inspector.vertex().right_child()) || mapping.count(inspector.vertex().right_child().negate())) {
      out << "(p) ";
    }
    PrintTo(inspector.vertex().right_child(), &out);
    out << std::endl;
    mapping.insert(std::make_pair(inspector.vertex(), true));

    ++inspector;
  }

  return out.str();
}

TEST(HashedReduce, StressTest) {
  const size_t REPEAT = 10000;
  CONSTEXPR_OR_CONST size_t RANK = 3;
  const size_t ENDOMORPHISMS_NUMBER = 20;

  typedef TVertexHashAlgorithms<hashers::SinglePowerHash, hashers::PermutationHash<Permutation16> > VertexHashAlgorithms;
  size_t seed = 0;
  while (++seed <= REPEAT) {
    UniformAutomorphismSLPGenerator<> generator(RANK, seed);
    auto image = EndomorphismSLP::composition(ENDOMORPHISMS_NUMBER, generator).image(1);

    Vertex reduced = VertexHashAlgorithms::reduce(image);

    std::vector<int> reduced_image;
    for (auto symbol : VertexWord(image)) {
      if (!reduced_image.empty() && symbol == -reduced_image.back()) {
        reduced_image.pop_back();
      } else {
        reduced_image.push_back(symbol);
      }
    }
    std::ostringstream reduced_image_string;
    std::copy(reduced_image.begin(), reduced_image.end(), std::ostream_iterator<int>(reduced_image_string, ""));
    auto correct_symbol = reduced_image.begin();
    for (auto symbol : VertexWord(reduced)) {
      ASSERT_EQ(*correct_symbol, symbol) << seed << std::endl
          << print_tree_preorder_single(image) << std::endl
          << print_tree_preorder(reduced) << std::endl
          << VertexWord(image) << std::endl
          << VertexWord(reduced) << std::endl
          << reduced_image_string.str() << std::endl;
      ++correct_symbol;
    }
  }
}

TEST(HashedReduceNarrow, StressTest) {
  const size_t REPEAT = 10000;
  CONSTEXPR_OR_CONST size_t RANK = 3;
  const size_t ENDOMORPHISMS_NUMBER = 20;

  typedef TVertexHashAlgorithms<hashers::SinglePowerHash, hashers::PermutationHash<Permutation16> > VertexHashAlgorithms;
  size_t seed = 0;
  while (++seed <= REPEAT) {
    UniformAutomorphismSLPGenerator<> generator(RANK, seed);
    auto image = EndomorphismSLP::composition(ENDOMORPHISMS_NUMBER, generator).image(1);

    Vertex reduced = VertexHashAlgorithms::reduce_narrow_slp(image);

    std::vector<int> reduced_image;
    for (auto symbol : VertexWord(image)) {
      if (!reduced_image.empty() && symbol == -reduced_image.back()) {
        reduced_image.pop_back();
      } else {
        reduced_image.push_back(symbol);
      }
    }
    std::ostringstream reduced_image_string;
    std::copy(reduced_image.begin(), reduced_image.end(), std::ostream_iterator<int>(reduced_image_string, ""));
    auto correct_symbol = reduced_image.begin();
    for (auto symbol : VertexWord(reduced)) {
      ASSERT_EQ(*correct_symbol, symbol) << seed << std::endl
          << print_tree_preorder_single(image) << std::endl
          << print_tree_preorder(reduced) << std::endl
          << VertexWord(image) << std::endl
          << VertexWord(reduced) << std::endl
          << reduced_image_string.str() << std::endl;
      ++correct_symbol;
    }
  }
}



TEST(HashedReduce, RemoveDuplicatesStressTest) {
  const size_t REPEAT = 1000;
  CONSTEXPR_OR_CONST size_t RANK = 3;
  const size_t ENDOMORPHISMS_NUMBER = 20;

  typedef TVertexHashAlgorithms<hashers::ImageLengthHash,
      hashers::SinglePowerHash,
      hashers::PermutationHash<Permutation16>
      > VertexHashAlgorithms;
  size_t seed = 0;
  size_t n_errors = 0;
  while (++seed <= REPEAT) {
    UniformAutomorphismSLPGenerator<> generator(RANK, seed);
    auto image = EndomorphismSLP::composition(ENDOMORPHISMS_NUMBER, generator).image(1);

    ASSERT_EQ(VertexHashAlgorithms::get_vertex_hash(image), VertexHashAlgorithms::get_vertex_hash(image));

    VertexHashAlgorithms::Cache cache;
    VertexHashAlgorithms::HashRepresentativesCache hash_cache;
    Vertex rd_mage = VertexHashAlgorithms::remove_duplicates(image, &cache, &hash_cache);

    ASSERT_EQ(image.length(), rd_mage.length());

    slp::MatchingTable mt;
    VertexWord image_word(image);
    VertexWord rd_word(rd_mage);
    if(!image_word.is_equal_to(rd_word, &mt)) {
      ++n_errors;
    }
  }
  ASSERT_EQ(0, n_errors);
}

} //anonymous namespace
} //slp
} //crag



