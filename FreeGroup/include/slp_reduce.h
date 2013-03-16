/*
 * slp_reduce.h
 *
 *  Created on: Mar 12, 2013
 *      Author: dpantele
 */

#pragma once
#ifndef CRAG_FREEGROUP_SLP_REDUCE_H_
#define CRAG_FREEGROUP_SLP_REDUCE_H_

#include <gmpxx.h>
typedef mpz_class LongInteger;

#include "slp_vertex.h"
#include "slp_mapper.h"
#include "slp_inspector.h"
#include "slp_common_prefix.h"

namespace std {
//(x)->_mp_d

template<>
struct hash<LongInteger> {
  private:
    constexpr static hash<mp_limb_t*> limb_hasher_ = hash<mp_limb_t*>();
  public:
    size_t operator()(const LongInteger& obj) const {
      return limb_hasher_(obj.get_mpz_t()->_mp_d);
    }
};

namespace tuple_hash_detail {

template <std::size_t I, typename Tuple>
class TupleHasher {
  private:
    constexpr static hash<typename std::tuple_element<I - 1, Tuple>::type> element_hasher_ = hash<typename std::tuple_element<I - 1, Tuple>::type>();
    constexpr static TupleHasher<I - 1, Tuple> prefix_hasher_ = TupleHasher<I - 1, Tuple>();
  public:
    size_t operator()(const Tuple& obj) const {
      size_t prefix_hash_value = prefix_hasher_(obj);
      //Taken from boost/functional/hash
      return element_hasher_(std::get<I - 1>(obj)) + 0x9e3779b9 + (prefix_hash_value << 6) + (prefix_hash_value >> 2);
    }
};

template <typename Tuple>
class TupleHasher<0, Tuple> {
  public:
    constexpr size_t operator()(const Tuple& obj) const {
      //Taken from boost/functional/hash
      return 0;
    }
};

template <std::size_t I, typename Tuple>
constexpr hash<typename std::tuple_element<I - 1, Tuple>::type> TupleHasher<I, Tuple>::element_hasher_;
template <std::size_t I, typename Tuple>
constexpr TupleHasher<I - 1, Tuple> TupleHasher<I, Tuple>::prefix_hasher_;
//template <typename Tuple>
//hash<typename std::tuple_element<0, Tuple>::type> TupleHasher<0, Tuple>::element_hasher_;

//template <class T>
//inline void hash_combine(std::size_t& seed, const T& v)
//{
//    hash<T> hasher;
//    seed ^= hasher(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
//}
//
//template <std::size_t I, typename T>
//inline typename std::enable_if<(I == std::tuple_size<T>::value),
//        void>::type
//hash_combine_tuple(std::size_t&, const T&)
//{ }
//
//template <std::size_t I, typename T>
//inline typename std::enable_if<(I < std::tuple_size<T>::value),
//        void>::type
//hash_combine_tuple(std::size_t& seed, const T& v)
//{
//    hash_combine(seed, std::get<I>(v));
//    hash_combine_tuple<I + 1>(seed, v);
//}
//
//template <typename T>
//inline std::size_t hash_tuple(T const& v)
//{
//    std::size_t seed = 0;
//    hash_combine_tuple<0>(seed, v);
//    return seed;
//}

}//namespace tuple_hash_detail

template<typename... T>
struct hash<tuple<T...>> {
  private:
    constexpr static tuple_hash_detail::TupleHasher<std::tuple_size<tuple<T...>>::value, tuple<T...>> tuple_haser_ = tuple_hash_detail::TupleHasher<std::tuple_size<tuple<T...>>::value, tuple<T...>>();
  public:
    size_t operator()(const tuple<T...>& obj) const {
      return tuple_haser_(obj);
    }
};

template<typename... T>
constexpr tuple_hash_detail::TupleHasher<std::tuple_size<tuple<T...>>::value, tuple<T...>> hash<tuple<T...>>::tuple_haser_;

}//namespace std

namespace crag {
namespace slp {

const Vertex& get_sub_slp(const Vertex& root, const LongInteger& begin, const LongInteger& end, std::unordered_map<std::tuple<Vertex, LongInteger, LongInteger>, Vertex>* cache);
inline Vertex get_sub_slp(const Vertex& root, const LongInteger& begin, const LongInteger& end) {
  std::unordered_map<std::tuple<Vertex, LongInteger, LongInteger>, Vertex> cache;
  return get_sub_slp(root, begin, end, &cache);
}

inline LongInteger get_cancellation_length(const Vertex& vertex, MatchingTable* matching_table) {
  return longest_common_prefix(vertex.left_child().negate(), vertex.right_child(), matching_table);
}

inline LongInteger get_cancellation_length(const Vertex& vertex) {
  MatchingTable temp;
  return get_cancellation_length(vertex, &temp);
}
inline Vertex reduce(const Vertex& vertex) {
  std::unordered_map<Vertex, Vertex> reduced_vertices;
  MatchingTable matching_table;
  map_vertices(vertex, &reduced_vertices, [&matching_table](const slp::Vertex& vertex, const std::unordered_map<Vertex, Vertex>& reduced_vertices) -> Vertex {
    if (vertex.height() <= 1) {
      return vertex;
    }
    const Vertex& left = reduced_vertices.find(vertex.left_child())->second;
    const Vertex& right = reduced_vertices.find(vertex.right_child())->second;
    if (!left && !right) {
      return Vertex();
    } else {
      NonterminalVertex result(left, right);
      LongInteger cancellation_length = get_cancellation_length(result, &matching_table);
      if (cancellation_length == 0) {
        if (left == vertex.left_child() && right == vertex.right_child()) {
          assert((vertex.height() > 1 && vertex.length() > 1) || (vertex.length() == 1 && vertex.height() == 1) || (vertex.length() == 0 && vertex.height() == 0));
          return vertex;
        } else if (!left) {
          return right;
        } else if (!right) {
          return left;
        } else {
          assert((result.height() > 1 && result.length() > 1) || (result.length() == 1 && result.height() == 1) || (result.length() == 0 && result.height() == 0));
          return result;
        }
      } else {
        Vertex reduced_left = get_sub_slp(left, 0, left.length() - cancellation_length);
        Vertex reduced_right = get_sub_slp(right, cancellation_length, right.length());
        if (!reduced_left && !reduced_right) {
          return Vertex();
        } else if (!reduced_left) {
          assert(reduced_right.height() >= 1);
          assert(reduced_right.height() != 1 || reduced_right.length() == 1);
          return reduced_right;
        } else if (!reduced_right) {
          assert(reduced_left.height() >= 1);
          assert(reduced_left.height() != 1 || reduced_left.length() == 1);
          return reduced_left;
        } else {
          auto result = NonterminalVertex(reduced_left, reduced_right);
          assert(result.length() > 1);
          assert(result.height() > 1);
          return result;
        }
      }
    }
  });
  return reduced_vertices[vertex];
}

}
}
#endif /* SLP_REDUCE_H_ */
