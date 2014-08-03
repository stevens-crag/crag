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

typedef crag::slp::TVertexHashAlgorithms<
    crag::slp::hashers::SinglePowerHash,
    crag::slp::hashers::PermutationHash<crag::Permutation16>
> VertexHashAlgorithms;

//const size_t LAST_REDUCTIONS_LIMIT = 7;
long int max_distance = 0;

std::map<long int, size_t> distances(long int distance = -1) {
  static std::map<long int, size_t> distances;
  if (distance != -1) {
    ++distances[distance];
  }

  return distances;
}

Vertex reduce(const Vertex& vertex,
              const std::unordered_map<Vertex, Vertex>& reduced_vertices,
              const LongInteger& left_siblings_length,
              std::map<LongInteger, size_t>* last_reductions,
              size_t this_reduction_id,
              VertexHashAlgorithms::Cache* hash_cache) {
  auto begin = std::chrono::system_clock::now();
  if (vertex.height() <= 1) {
    return vertex;
  }
  auto reversed = reduced_vertices.find(vertex.negate());
  if (reversed != reduced_vertices.end()) {
    return reversed->second.negate();
  }
  const Vertex& left = reduced_vertices.find(vertex.left_child())->second;
  const Vertex& right = reduced_vertices.find(vertex.right_child())->second;
  if (!left && !right) {
    return Vertex();
  } else {
    NonterminalVertex result(left, right);

    bool success = false;
    auto first = left.negate();
    auto second = right;

    if (first.length() > second.length()) {
      std::swap(first, first);
    }

    LongInteger cancellation_length = 0;

    auto first_prefix = first;
    auto second_prefix = second;

    while (first_prefix.height() > 0 && second_prefix.height() > 0) {
      if (first_prefix == second_prefix) {
        cancellation_length = first_prefix.length();
        if (get_sub_slp(first, first_prefix.length(), first_prefix.length() + 1) != get_sub_slp(second, first_prefix.length(), first_prefix.length() + 1)) {
          success = true;
        }
        break;
      }
      if (first_prefix.length() >= second_prefix.length()) {
        first_prefix = first_prefix.left_child();
      } else {
        second_prefix = second_prefix.left_child();
      }
    }

    if (first_prefix.height() == 0 || second_prefix.height() == 0) {
      success = true;
    }

    if (!success) {
      for (auto& reduction : *last_reductions) {
        if (reduction.first <= cancellation_length) {
          continue;
        }

        auto first_hash = VertexHashAlgorithms::get_subvertex_hash(first, 0, reduction.first, hash_cache);

        if (first_hash == VertexHashAlgorithms::get_subvertex_hash(second, 0, reduction.first, hash_cache)) {
          cancellation_length = reduction.first;
          if (crag::slp::get_sub_slp(first, cancellation_length, cancellation_length + 1) != crag::slp::get_sub_slp(second, cancellation_length, cancellation_length + 1)) {
            auto last_id = reduction.second;
            auto current_distance = std::count_if(last_reductions->begin(), last_reductions->end(), [last_id](const std::pair<const LongInteger, size_t>& length) {
              return length.second > last_id;
            });

            max_distance = std::max(max_distance, current_distance);
            distances(current_distance);

            reduction.second = this_reduction_id;
            success = true;
            break;
          }
        } else {
          break;
        }
      }
    }

    if (!success) {
      cancellation_length = VertexHashAlgorithms::get_longest_common_prefix(first, second, cancellation_length, -1, hash_cache);
//      LongInteger cancellation_length = VertexHashAlgorithms::get_cancellation_length(result, hash_cache);
      (*last_reductions)[cancellation_length] = this_reduction_id;
    }

    if (cancellation_length == 0) {
      if (left == vertex.left_child() && right == vertex.right_child()) {
        assert(
            (vertex.height() > 1 && vertex.length() > 1) || (vertex.length() == 1 && vertex.height() == 1) || (vertex.length() == 0 && vertex.height() == 0));
        auto end = std::chrono::system_clock::now();
        return vertex;
      } else if (!left) {
        auto end = std::chrono::system_clock::now();
        return right;
      } else if (!right) {
        auto end = std::chrono::system_clock::now();
        return left;
      } else {
        auto end = std::chrono::system_clock::now();
        assert(
            (result.height() > 1 && result.length() > 1) || (result.length() == 1 && result.height() == 1) || (result.length() == 0 && result.height() == 0));
        return result;
      }
    } else {

      Vertex reduced_left = get_sub_slp(left, 0,
                                        left.length() - cancellation_length);
      Vertex reduced_right = get_sub_slp(right, cancellation_length,
                                         right.length());

      if (!reduced_left && !reduced_right) {
        auto end = std::chrono::system_clock::now();
        return Vertex();
      } else if (!reduced_left) {
        assert(reduced_right.height() >= 1);
        assert(reduced_right.height() != 1 || reduced_right.length() == 1);
        auto end = std::chrono::system_clock::now();
        return reduced_right;
      } else if (!reduced_right) {
        assert(reduced_left.height() >= 1);
        assert(reduced_left.height() != 1 || reduced_left.length() == 1);
        auto end = std::chrono::system_clock::now();
        return reduced_left;
      } else {
        auto result = NonterminalVertex(reduced_left, reduced_right);
        auto end = std::chrono::system_clock::now();
        assert(result.length() > 1);
        assert(result.height() > 1);
        return result;
      }
    }
  }
}

int main() {
  gmp_pool_setup();
  size_t seed = 112233;

  typedef crag::slp::TVertexHashAlgorithms<
      crag::slp::hashers::SinglePowerHash,
      crag::slp::hashers::PermutationHash<crag::Permutation16>
  > VertexHashAlgorithms;

  size_t endomorphisms_number = 10;
  while (endomorphisms_number < 10000) {
    for (auto rank : {3, 6, 10, 20, 40, 100}) {
      UniformAutomorphismSLPGenerator<> generator(rank, seed + endomorphisms_number + rank);
      auto slp = EndomorphismSLP::composition(endomorphisms_number, generator).image(1);
      max_distance = 0;

      std::unordered_map<Vertex, Vertex> reduced_vertices;
      VertexHashAlgorithms::Cache hash_cache;

      auto acceptor = [&reduced_vertices] (const crag::slp::inspector::InspectorTask& task) {
        return reduced_vertices.count(task.vertex) == 0;//do not accept if vertex is mapped already
      };

      crag::slp::Inspector<crag::slp::inspector::Postorder, decltype(acceptor)> inspector(slp, acceptor);

      size_t reduce_count = 0;
      std::map<LongInteger, size_t> last_reductions;

      auto begin = std::chrono::system_clock::now();

      while (!inspector.stopped()) {
        auto img = reduce(inspector.vertex(), reduced_vertices, inspector.vertex_left_siblings_length(), &last_reductions, ++reduce_count, &hash_cache);
        auto new_entry = std::make_pair(inspector.vertex(), img);

        reduced_vertices.insert(new_entry);
        inspector.next();
      }
      int total = 0;
      std::cout << max_distance << std::endl;
      for (auto distance : distances()) {
        std::cout << distance.first << ' ' << distance.second << ' ' << (total += distance.second) << std::endl;
      }

      std::cout << "End#: " << endomorphisms_number << ", rank: " << rank << ", max_distance: " << max_distance << std::endl;
    }
    endomorphisms_number *= 2;
  }
  return 0;
}


