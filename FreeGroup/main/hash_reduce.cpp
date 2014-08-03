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
//#include "boost/functional/hash/extensions.hpp"

using crag::slp::Vertex;
using crag::slp::VertexWord;
using crag::slp::PreorderInspector;
using crag::slp::TerminalVertex;
using crag::slp::NonterminalVertex;
using crag::slp::MatchingTable;
using crag::UniformAutomorphismSLPGenerator;
using crag::EndomorphismSLP;

//namespace crag { namespace slp {
//std::string print_tree_preorder_single(const Vertex& vertex) {
//  std::ostringstream out;
//  std::unordered_map<slp::Vertex, bool> mapping;
//  mapper::SkipMappedVerticesAcceptor<std::unordered_map<slp::Vertex, bool>> acceptor(mapping);
//  Inspector<inspector::Preorder, mapper::SkipMappedVerticesAcceptor<std::unordered_map<slp::Vertex, bool>>> inspector(vertex, acceptor);
//  while (!inspector.stopped()) {
//    PrintTo(inspector.vertex(), &out);
//    out << " << ";
//    if (mapping.count(inspector.vertex().left_child()) || mapping.count(inspector.vertex().left_child().negate())) {
//      out << "(p) ";
//    }
//    PrintTo(inspector.vertex().left_child(), &out);
//    out << " >> ";
//    if (mapping.count(inspector.vertex().right_child()) || mapping.count(inspector.vertex().right_child().negate())) {
//      out << "(p) ";
//    }
//    PrintTo(inspector.vertex().right_child(), &out);
//    out << std::endl;
//    mapping.insert(std::make_pair(inspector.vertex(), true));
//
//    ++inspector;
//  }
//
//  return out.str();
//}
//
//}}

int main() {
  gmp_pool_setup();
  int REPEAT = 1;
  CONSTEXPR_OR_CONST size_t RANK = 6;
  CONSTEXPR_OR_CONST size_t ENDOMORPHISMS_NUMBER = 2200;
  size_t seed = 112233;
  UniformAutomorphismSLPGenerator<> generator(RANK, seed);
  auto begin = std::chrono::system_clock::now();
  int count = REPEAT;

  LongInteger wrong_answers;
  LongInteger total_calculations;
  typedef crag::slp::TVertexHashAlgorithms<
      //crag::slp::hashers::PowerCountHash<RANK>,
      crag::slp::hashers::PermutationHash<crag::Permutation16>
  > VertexHashAlgorithms;

  typedef crag::slp::TVertexHashAlgorithms<
      crag::slp::hashers::SinglePowerHash,
      crag::slp::hashers::PermutationHash<crag::Permutation16>
  > WeakVertexHashAlgorithms;

  auto strong_duration = begin - begin;
  auto weak_duration = begin - begin;

  size_t reduce_different = 0;

  while (--count >= 0) {
    auto image = EndomorphismSLP::composition(ENDOMORPHISMS_NUMBER, generator).image(1);
    crag::slp::MatchingTable matching_table;
    std::unordered_map<Vertex, Vertex> reduced_vertices;
    WeakVertexHashAlgorithms::Cache vertex_hashes;
    auto reduce_start = std::chrono::system_clock::now();
    //Vertex reduced = VertexHashAlgorithms::reduce(image);
    strong_duration += std::chrono::system_clock::now() - reduce_start;
    reduce_start = std::chrono::system_clock::now();
    Vertex weak_reduced = WeakVertexHashAlgorithms::reduce(image, &vertex_hashes, &reduced_vertices);
    weak_duration += std::chrono::system_clock::now() - reduce_start;
    //if (reduced.length() != weak_reduced.length()) {
      //++reduce_different;
    //}
    std::unordered_map<WeakVertexHashAlgorithms::VertexHash, size_t> hash_count;
    std::unordered_set<Vertex> examined_vertices;
    auto acceptor = [&examined_vertices] (const crag::slp::inspector::InspectorTask& task) -> bool {
      return !examined_vertices.count(task.vertex);//do not accept if vertex is mapped already
    };

    crag::slp::Inspector<crag::slp::inspector::Preorder, decltype(acceptor)> original_inspector(image, acceptor);

    size_t vertex_count = 0;
    while (!original_inspector.stopped()) {
      ++hash_count[WeakVertexHashAlgorithms::get_subvertex_hash(original_inspector.vertex(), 0, original_inspector.vertex().length(), &vertex_hashes)];
      examined_vertices.insert(original_inspector.vertex());
      original_inspector.next();
      ++vertex_count;
    }
    std::cout << "Vertices: " << vertex_count << std::endl;

    std::map<size_t, size_t> count_count;

    for (auto count : hash_count) {
      ++count_count[count.second];
    }

    for (auto count : count_count) {
      std::cout << count.first << ": " << count.second << std::endl;
    }

    hash_count.clear();
    examined_vertices.clear();

    crag::slp::Inspector<crag::slp::inspector::Preorder, decltype(acceptor)> inspector(weak_reduced, acceptor);

    vertex_count = 0;
    while (!inspector.stopped()) {
      ++hash_count[WeakVertexHashAlgorithms::get_subvertex_hash(inspector.vertex(), 0, inspector.vertex().length(), &vertex_hashes)];
      examined_vertices.insert(inspector.vertex());
      inspector.next();
      ++vertex_count;
    }
    std::cout << "Vertices: " << vertex_count << std::endl;

    count_count.clear();

    for (auto count : hash_count) {
      ++count_count[count.second];
    }

    for (auto count : count_count) {
      std::cout << count.first << ": " << count.second << std::endl;
    }

    auto end = std::chrono::system_clock::now();
    std::cout << "Duration: "<< std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << std::endl;
  }

  auto end = std::chrono::system_clock::now();
  std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()/REPEAT << std::endl;
  std::cout << reduce_different << " of " << REPEAT << std::endl;
  std::cout << "Strong: " << std::chrono::duration_cast<std::chrono::milliseconds>(strong_duration).count()/REPEAT << std::endl;
  std::cout << "Weak: " << std::chrono::duration_cast<std::chrono::milliseconds>(weak_duration).count()/REPEAT << std::endl;

}


