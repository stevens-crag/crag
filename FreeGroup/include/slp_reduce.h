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

inline LongInteger get_cancellation_length(const Vertex& vertex) {
  return longest_common_prefix(vertex.left_child(), vertex.right_child().negate());
}

inline Vertex reduce(const Vertex& vertex) {
  std::unordered_map<Vertex, Vertex> reduced_vertices;
  map_vertices(vertex, &reduced_vertices, [](const slp::Vertex& vertex, std::unordered_map<Vertex, Vertex>* reduced_vertices) -> Vertex {
    const Vertex& left = (*reduced_vertices)[vertex.left_child()];
    const Vertex& right = (*reduced_vertices)[vertex.right_child()];
    if (!left && !right) {
      return Vertex();
    } else {
      NonterminalVertex result(left, right);
      LongInteger cancellation_length = get_cancellation_length(result);
      if (cancellation_length == 0) {
        if (left == vertex.left_child() && right == vertex.right_child()) {
          return vertex;
        } else {
          return result;
        }
      } else {
        Vertex reduced_left = get_sub_slp(left, 0, left.length() - cancellation_length);
        Vertex reduced_right = get_sub_slp(right, cancellation_length, right.length());
        if (!reduced_left && !reduced_right) {
          return Vertex();
        } else if (!reduced_left) {
          return reduced_right;
        } else if (!reduced_right) {
          return reduced_left;
        } else {
          return NonterminalVertex(reduced_left, reduced_right);
        }
      }
    }
  });
  return reduced_vertices[vertex];
}

}
}
#endif /* SLP_REDUCE_H_ */
