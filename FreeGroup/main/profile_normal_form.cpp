/*
 * profile_matching_new.cpp
 *
 *  Created on: Apr 2, 2013
 *      Author: dpantele
 */

#include <memory>
#include <functional>
#include <utility>
#include <vector>
#include <tuple>
#include <algorithm>
#include <chrono>
#include <array>

#include "slp.h"
#include "EndomorphismSLP.h"
#include "Permutation.h"
#include "gmp_boost_pool_allocator.h"
#include "boost/functional/hash/extensions.hpp"

using crag::slp::Vertex;
using crag::slp::VertexWord;
using crag::slp::PreorderInspector;
using crag::slp::TerminalVertexTemplate;
typedef crag::slp::TerminalVertexTemplate<int> TerminalVertex;
using crag::slp::NonterminalVertex;
using crag::slp::MatchingTable;
using crag::UniformAutomorphismSLPGenerator;
using crag::EndomorphismSLP;

int main() {
  gmp_pool_setup();
  int REPEAT = 5;
  constexpr size_t RANK = 6;
  constexpr size_t ENDOMORPHISMS_NUMBER = 2000;
  size_t seed = 112233;
  UniformAutomorphismSLPGenerator<int> generator(RANK, seed);
  auto begin = std::chrono::system_clock::now();
  int count = REPEAT;

//  typedef crag::slp::TVertexHashAlgorithms<
//      crag::slp::hashers::SinglePowerHash,
//      crag::slp::hashers::PermutationHash<crag::Permutation16>
//  > VertexHashAlgorithms;

  auto normal_form_duration = begin - begin;

  while (--count >= 0) {
    auto image = EndomorphismSLP<int>::composition(ENDOMORPHISMS_NUMBER, generator).image(1);

//    std::unordered_map<Vertex, Vertex> reduced_vertices;
//    VertexHashAlgorithms::Cache vertex_hashes;
//    Vertex reduced = VertexHashAlgorithms::reduce(image, &vertex_hashes, &reduced_vertices);

    auto normal_form_start = std::chrono::system_clock::now();
    Vertex normal_form = crag::slp::recompression::normal_form(image);
    normal_form_duration += std::chrono::system_clock::now() - normal_form_start;
    auto end = std::chrono::system_clock::now();
    std::cout << "Duration: "<< std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << std::endl;
  }

  auto end = std::chrono::system_clock::now();
  std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()/REPEAT << std::endl;
  std::cout << "Normal form: " << std::chrono::duration_cast<std::chrono::milliseconds>(normal_form_duration).count()/REPEAT << std::endl;
}


