/*
 * EndomorphsimSLP.h
 *
 *  Created on: Jan 29, 2013
 *      Author: pmorar
 */

#ifndef CRAG_FREE_GROUPS_ENDOMORPHISM_SLP_H_
#define CRAG_FREE_GROUPS_ENDOMORPHISM_SLP_H_

#include <map>
#include <unordered_map>
#include <random>
#include <assert.h>
#include "slp.h"

namespace crag {

/**
 * Represents a free group endomorphism using straight-line programs.
 * @tparam TerminalSymbol terminal symbols class
 */
template<typename TerminalSymbol = int>
class EndomorphismSLP {
public:

  typedef typename std::map<TerminalSymbol, slp::Vertex>::size_type size_type;

  //use default copy/move constructors/assignments

  //! Returns the identity automorphism
  static EndomorphismSLP identity() {
    return EndomorphismSLP();
  }

  //! Returns the automorphism inverting the specified terminal symbbol.
  /**
   * @param symbol inverted terminal symbol, must be > 0
   */
  static EndomorphismSLP inverter(const TerminalSymbol& symbol) {
    assert(is_positive_terminal_symbol(symbol));
    return EndomorphismSLP(symbol);
  }

  //! Returns the automorphism mapping #symbol to the product #symbol #right_multiplier.
  /**
   * @param symbol mapped terminal symbol, must be > 0
   * @param	right_multiplier	right multiplier
   */
  static EndomorphismSLP right_multiplier(const TerminalSymbol& symbol, const TerminalSymbol& right_multiplier) {
    assert(is_positive_terminal_symbol(symbol));
    assert(symbol != right_multiplier);
    EndomorphismSLP tmp;
    tmp.images_.insert(std::make_pair(symbol,
      slp::NonterminalVertex(TerminalVertex(symbol), TerminalVertex(right_multiplier))));
    return tmp;
  }

  //! Returns the automorphism mapping #symbol to the product #left_multiplier #symbol.
  /**
   * @param  left_multiplier  left multiplier
   * @param symbol mapped terminal symbol, must be > 0
   */
  static EndomorphismSLP left_multiplier(const TerminalSymbol& left_multiplier, const TerminalSymbol& symbol) {
    assert(is_positive_terminal_symbol(symbol));
    assert(left_multiplier != symbol);
    EndomorphismSLP tmp;
    tmp.images_.insert(std::make_pair(symbol,
      slp::NonterminalVertex(TerminalVertex(left_multiplier), TerminalVertex(symbol))));
    return tmp;
  }

  //! Returns the composition of endomorphisms specified by the range
  template<typename Iterator>
  static EndomorphismSLP composition(Iterator begin, Iterator end) {
    return identity().compose_with(begin, end);
  }

  //! Returns the composition of #num endomorphisms produced by #generator
  /**
   * @param num      the number of endomorphisms to compose
   * @param generator endomorphisms generator (using operator())
   */
  template<typename Generator>
  static EndomorphismSLP composition(unsigned int num, Generator& generator) {
    return identity().compose_with(num, generator);
  }

  //! Compose with the given endomorphism.
  EndomorphismSLP& operator*=(const EndomorphismSLP& a);

  //! Compose with endomorphisms specified by the range.
  template<typename Iterator>
  EndomorphismSLP& compose_with(Iterator begin, Iterator end) {
    for(;begin != end; ++begin)
      (*this) *= *begin;
    return *this;
  }

  //! Compose with #num automorphism produced by #generator
  /**
   * @param num the number of endomorphisms to compose with
   * @param generator endomorphisms generator (using operator())
   */
  template<typename Generator>
  EndomorphismSLP& compose_with(unsigned int num, Generator& generator) {
    for (unsigned int i = 0; i < num; ++i)
      (*this) *= generator();
    return *this;
  }


  //! Returns the image of the terminal.
  slp::VertexWord<TerminalSymbol> image(const TerminalSymbol& t) const {
    bool is_positive = is_positive_terminal_symbol(t);
    return slp::VertexWord<TerminalSymbol>(is_positive ? slp(t) : slp(-t).negate());
  }

  //! Returns the root of the straight-line program representing the positive terminal symbol.
  /**
   * @param t positive terminal symbol
   */
  slp::Vertex slp(const TerminalSymbol& t) const {
    assert(is_positive_terminal_symbol(t));
    auto result = images_.find(t);
    if (result == images_.end()) //if it is not in images_, then it is the identity map.
      return TerminalVertex(t);
    else
      return result->second;
  }

  //! Returns the maximal terminal symbol with non-identity image.
  TerminalSymbol max_non_trivial_image_symbol() const {
    if (images_.empty())
      return TerminalSymbol();
    return images_.rbegin()->first;
  }

  size_type non_trivial_images_num() const {
    return images_.size();
  }

private:
  typedef typename slp::TerminalVertexTemplate<TerminalSymbol> TerminalVertex;

  //! Checks whether the symbol is not inverse.
  static bool is_positive_terminal_symbol(const TerminalSymbol& symbol) {
    static const TerminalSymbol null_symbol = TerminalSymbol();
    return symbol > null_symbol;
  }

  //! Helper method that copies #vertex such that the children of the new vertex are images of the children of #vertex
  /**
   * It assumes that if #vertex is non-terminal, then #images contains its children images
   */
  slp::Vertex map_vertex(const slp::Vertex& vertex, const std::unordered_map<slp::Vertex, slp::Vertex>& images) const;

  //! The default constructor.
  EndomorphismSLP() {}

  EndomorphismSLP(const TerminalSymbol& inverted) {
    images_.insert(std::make_pair(inverted,
      TerminalVertex(inverted).negate()));
  }

  //! The representation of images of terminal symbols by straight-line programs
  /**
   * If there is no entry for a given terminal symbol, then its image is the terminal itself.
   * Also we use {@code std::map} instead of {@code std::unordered_map} because we would like
   * to iterate over the keys in the specific order defined by operator < for TerminalSymbol.
   */
  std::map<TerminalSymbol, slp::Vertex> images_;

  class MapperBinder {
  public:
    MapperBinder(const EndomorphismSLP& endomorphism)
      : endomorphism_(endomorphism){
    }

    slp::Vertex operator()(const slp::Vertex& vertex, const std::unordered_map<slp::Vertex, slp::Vertex>& images) const {
      return endomorphism_.map_vertex(vertex, images);
    }

  private:
    const EndomorphismSLP& endomorphism_;
  };
};

//! Compose given endomorphisms.
/**
 * Returns const to behave consistently with built-in types.
 */
template<typename TerminalSymbol>
const EndomorphismSLP<TerminalSymbol> operator*(const EndomorphismSLP<TerminalSymbol>& e1, const EndomorphismSLP<TerminalSymbol>& e2) {
  return EndomorphismSLP<TerminalSymbol>(e1) *= e2;
}



//! Automorphisms generator
/**
 * @tparam TerminalSymbol          terminal symbol representation. We suppose it has a constructor taking index
 * @tparam RandomEngine            engine generating uniformly random non-negative numbers. See std::random library documentation.
 * @tparam TerminalSymbolIndexType terminal symbol index type
 */
template <typename TerminalSymbol = int,
  typename RandomEngine = std::default_random_engine>
class UniformAutomorphismSLPGenerator {
public:
  typedef int index_type;

  //! Constructs a generator of automorphisms of the free group of the given rank.
  /**
   * @param rank free group rank > 0
   */
  UniformAutomorphismSLPGenerator(index_type rank)
    : UniformAutomorphismSLPGenerator(rank, std::make_shared<RandomEngine>(), nullptr)
  {}

  //! Constructs a generator of automorphisms of the free group of the given rank.
  /**
   * @param rank free group rank > 0
   * @param seed random engine seed for creation of a new one
   */
  explicit UniformAutomorphismSLPGenerator(index_type rank, typename RandomEngine::result_type seed)
    : UniformAutomorphismSLPGenerator(rank, ::std::make_shared<RandomEngine>(seed), nullptr)
  {}

  //! Constructs a generator of automorphisms of the free group of the given rank.
  /**
   * @param rank free group rank > 0
   * @param random_engine random engine
   */
  explicit UniformAutomorphismSLPGenerator(index_type rank, RandomEngine* random_engine)
    : UniformAutomorphismSLPGenerator(rank, nullptr, random_engine)
  {}

  //! Generates a random automorphism
  EndomorphismSLP<TerminalSymbol> operator()() {
    index_type r_val = random_distr_(*random_engine_);
    if (r_val < MULTIPLIERS_COUNT) {
      const bool right_multiplier = (r_val % 2) == 0;
      r_val >>= 1;
      const int mapped_symbol_index = 1 + ( r_val % RANK );
    const TerminalSymbol mapped_symbol(mapped_symbol_index);
    const int multiplier_index = 1 + ( r_val / RANK );
    const TerminalSymbol multiplier(multiplier_index < mapped_symbol_index ? multiplier_index : multiplier_index + 1);

    if (right_multiplier)
      return EndomorphismSLP<TerminalSymbol>::right_multiplier(mapped_symbol, multiplier);
    else
      return EndomorphismSLP<TerminalSymbol>::left_multiplier(multiplier, mapped_symbol);

  } else { // r_val >= MULTIPLIERS_COUNT
    return EndomorphismSLP<TerminalSymbol>::inverter(TerminalSymbol(1 + r_val - MULTIPLIERS_COUNT));
  }
  }


private:
  const index_type RANK;
  const index_type RIGHT_MULTIPLIERS_COUNT;
  const index_type LEFT_MULTIPLIERS_COUNT = RIGHT_MULTIPLIERS_COUNT;
  const index_type MULTIPLIERS_COUNT = RIGHT_MULTIPLIERS_COUNT + LEFT_MULTIPLIERS_COUNT;
  const index_type INVERTERS_COUNT;
  const index_type COUNT = MULTIPLIERS_COUNT + INVERTERS_COUNT;

  UniformAutomorphismSLPGenerator(index_type rank, const ::std::shared_ptr<RandomEngine>& random_engine_ptr, RandomEngine* random_engine)
      : RANK(rank),
        RIGHT_MULTIPLIERS_COUNT(rank * (rank - 1)),
        INVERTERS_COUNT(rank),
        random_engine_ptr_(random_engine_ptr),
        random_engine_(random_engine_ptr ? random_engine_ptr.get() : random_engine),
        random_distr_(0, COUNT - 1) {
      assert(rank > 0);
    }

  ::std::shared_ptr<RandomEngine> random_engine_ptr_;
  RandomEngine* random_engine_;
  std::uniform_int_distribution<index_type> random_distr_;
};


template <typename TerminalSymbol>
EndomorphismSLP<TerminalSymbol>& EndomorphismSLP<TerminalSymbol>::operator*=(const EndomorphismSLP<TerminalSymbol>& a) {
  std::unordered_map<slp::Vertex, slp::Vertex> new_vertices;//a's vertices to new vertices correspondence

  MapperBinder binder(*this);
  for (auto root_entry: a.images_) {//mapping vertices of #a to new ones
    slp::map_vertices(root_entry.second, &new_vertices, binder);
  }

  //replacing roots
  std::map<TerminalSymbol, slp::Vertex> new_images;
  for (auto root_entry: a.images_) {
    auto new_root = new_vertices.find(root_entry.second)->second;
      new_images.insert(std::make_pair(root_entry.first, new_root));
  }
  //adding images that were not inserted
  for (auto root_entry: images_) {
    if (new_images.find(root_entry.first) == new_images.end())//it was not mapped by a
      new_images.insert(root_entry);
  }

  swap(images_, new_images);
  return *this;
}

template<typename TerminalSymbol>
slp::Vertex EndomorphismSLP<TerminalSymbol>::map_vertex(const slp::Vertex& vertex, const std::unordered_map<slp::Vertex, slp::Vertex>& images) const {
  if (!vertex)
    return vertex;//Mapping null vertex

  if (vertex.height() == 1) {//the vertex is terminal so map it to our corresponding root
    const TerminalSymbol& symbol = TerminalVertex(vertex).terminal_symbol();
    bool is_positive = is_positive_terminal_symbol(symbol);
    auto positive_symbol = is_positive ? symbol : -symbol;
    slp::Vertex v = slp(positive_symbol);
    if (TerminalVertex(positive_symbol) == v)//if id map
      return vertex;
    return is_positive ? v : v.negate();
  } else {//for a nonterminal we already processed its children because postorder inspector
    auto left_val = images.find(vertex.left_child());
    auto right_val = images.find(vertex.right_child());
    assert(left_val != images.end() && right_val != images.end());
    auto left = left_val->second;
    auto right = right_val->second;
    if (left == vertex.left_child() && right == vertex.right_child()) //if children were not copied, then we should not copy vertex
      return vertex;
    return slp::NonterminalVertex(left, right);
  }
}

} /* namespace crag */
#endif /* CRAG_FREE_GROUPS_ENDOMORPHISM_SLP_H_ */
