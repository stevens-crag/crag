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

Vertex logging_get_sub_slp(const Vertex& root, const LongInteger& begin, const LongInteger& end) {
  if (begin >= root.length() || end < 0 || end <= begin) {
    static Vertex Null;
    return Null;
  }
  if (root.height() == 1) {
    assert(root.length() == 1);
    return root;
  } else if (begin <= 0 && end >= root.length()) {
    assert(root.height() != 1 || root.length() == 1);
    return root;
  } else {
    if (root.split_point() >= end) {
      std::cout << root.right_child().vertex_id() << " (r" << (root.right_child().height() > 1 ? "" : "t") << ") * ";
      return logging_get_sub_slp(root.left_child(), begin, end);
    } else if (root.split_point() <= begin) {
      std::cout << root.left_child().vertex_id() << " (l" << (root.left_child().height() > 1 ? "" : "t") << ") * ";
      return logging_get_sub_slp(root.right_child(), begin - root.split_point(), end - root.split_point());
    } else {
      return NonterminalVertex(
          logging_get_sub_slp(root.left_child(), begin, root.split_point()),
          logging_get_sub_slp(root.right_child(), 0, end - root.split_point())
      );
    }
  }
}

const size_t LAST_REDUCTIONS_LIMIT = 7;
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
      if (last_reductions->size() >= LAST_REDUCTIONS_LIMIT) {
        auto oldest_entry = min_element(last_reductions->begin(), last_reductions->end(),
            [](const std::pair<const LongInteger, size_t>& first, const std::pair<const LongInteger, size_t>& second) {
              return first.second < second.second;
            }
        );

        last_reductions->erase(oldest_entry);
      }

      (*last_reductions)[cancellation_length] = this_reduction_id;
    }

    std::cout << vertex.vertex_id() << ',' << vertex.length() << ',' << vertex.height() << ',' << vertex.left_child().vertex_id() << ',' << vertex.right_child().vertex_id() << ',' << 2*cancellation_length << ',' << left_siblings_length << ',' << success << ',';
    if (cancellation_length == 0) {
      if (left == vertex.left_child() && right == vertex.right_child()) {
        assert(
            (vertex.height() > 1 && vertex.length() > 1) || (vertex.length() == 1 && vertex.height() == 1) || (vertex.length() == 0 && vertex.height() == 0));
        auto end = std::chrono::system_clock::now();
        std::cout << ',' << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
        std::cout << "=\n";
        return vertex;
      } else if (!left) {
        auto end = std::chrono::system_clock::now();
        std::cout << ',' << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
        std::cout << ">\n";
        return right;
      } else if (!right) {
        auto end = std::chrono::system_clock::now();
        std::cout << ',' << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
        std::cout << "<\n";
        return left;
      } else {
        auto end = std::chrono::system_clock::now();
        std::cout << ',' << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
        std::cout << "!\n";
        assert(
            (result.height() > 1 && result.length() > 1) || (result.length() == 1 && result.height() == 1) || (result.length() == 0 && result.height() == 0));
        return result;
      }
    } else {

      Vertex reduced_left = logging_get_sub_slp(left, 0,
                                        left.length() - cancellation_length);
      std::cout << ',';
      Vertex reduced_right = logging_get_sub_slp(right, cancellation_length,
                                         right.length());

      if (!reduced_left && !reduced_right) {
        auto end = std::chrono::system_clock::now();
        std::cout << ',' << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
        std::cout << ", 0, e\n";
        return Vertex();
      } else if (!reduced_left) {
        assert(reduced_right.height() >= 1);
        assert(reduced_right.height() != 1 || reduced_right.length() == 1);
        auto end = std::chrono::system_clock::now();
        std::cout << ',' << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
        std::cout << ',' << reduced_right.vertex_id() << ", r\n";
        return reduced_right;
      } else if (!reduced_right) {
        assert(reduced_left.height() >= 1);
        assert(reduced_left.height() != 1 || reduced_left.length() == 1);
        auto end = std::chrono::system_clock::now();
        std::cout << ',' << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
        std::cout << ',' << reduced_left.vertex_id()<< ", l\n";
        return reduced_left;
      } else {
        auto result = NonterminalVertex(reduced_left, reduced_right);
        auto end = std::chrono::system_clock::now();
        std::cout << ',' << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
//        auto temp_result_joined = NonterminalVertex(result, result);
//        auto conjugate_length = VertexHashAlgorithms::get_cancellation_length(temp_result_joined, hash_cache);
//        std::cout << ',' << result.vertex_id() << ',' << conjugate_length << ", n\n";
        std::cout << ',' << result.vertex_id() << ", n\n";
        assert(result.length() > 1);
        assert(result.height() > 1);
        return result;
      }
    }
  }
}

int main() {
  gmp_pool_setup();
  CONSTEXPR_OR_CONST size_t RANK = 6;
  CONSTEXPR_OR_CONST size_t ENDOMORPHISMS_NUMBER = 2200;
  size_t seed = 112233;
  UniformAutomorphismSLPGenerator<> generator(RANK, seed);

  typedef crag::slp::TVertexHashAlgorithms<
      crag::slp::hashers::SinglePowerHash,
      crag::slp::hashers::PermutationHash<crag::Permutation16>
  > VertexHashAlgorithms;

  auto slp = EndomorphismSLP::composition(ENDOMORPHISMS_NUMBER, generator).image(1);

  std::unordered_map<Vertex, Vertex> reduced_vertices;
  VertexHashAlgorithms::Cache hash_cache;

  auto acceptor = [&reduced_vertices] (const crag::slp::inspector::InspectorTask& task) {
    return reduced_vertices.count(task.vertex) == 0;//do not accept if vertex is mapped already
  };

  crag::slp::Inspector<crag::slp::inspector::Postorder, decltype(acceptor)> inspector(slp, acceptor);


  size_t reduce_count;
  std::map<LongInteger, size_t> last_reductions;

  auto begin = std::chrono::system_clock::now();

  while (!inspector.stopped()) {
    auto& vertex = inspector.vertex();

    auto img = reduce(inspector.vertex(), reduced_vertices, inspector.vertex_left_siblings_length(), &last_reductions, ++reduce_count, &hash_cache);
    auto new_entry = std::make_pair(inspector.vertex(), img);

    reduced_vertices.insert(new_entry);
    inspector.next();
  }
  auto end = std::chrono::system_clock::now();

  for (auto& item : last_reductions) {
    std::cout << item.first << ',' << item.second << std::endl;
  }

  std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << std::endl;

  return 0;
}


