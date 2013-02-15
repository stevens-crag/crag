/*
 * EndomorphsimSLP.h
 *
 *  Created on: Jan 29, 2013
 *      Author: pmorar
 */

#ifndef CRAG_FREE_GROUPS_ENDOMORPHISM_SLP_H_
#define CRAG_FREE_GROUPS_ENDOMORPHISM_SLP_H_

#include <map>
#include <random>
#include <assert.h>
#include "SLPSet.h"

namespace crag {

/**
 * Represents a free group endomorphism using straight-line programs.
 * @tparam TerminalSymbol terminal symbols class
 */
template<typename TerminalSymbol = int>
class EndomorphismSLP {
public:

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
	  assert(symbol > 0);
		return EndomorphismSLP(symbol);
	}

	//! Returns the automorphism mapping #symbol to the product #symbol #right_multiplier.
	/**
	 * @param symbol mapped terminal symbol, must be > 0
	 * @param	right_multiplier	right multiplier
	 */
	static EndomorphismSLP right_multiplier(const TerminalSymbol& symbol, const TerminalSymbol& right_multiplier) {
	  assert(symbol > 0);
	  EndomorphismSLP tmp;
	  auto image_vertex = SLPVertex::concatenate(
      SLPVertex::terminal_vertex(symbol),
      SLPVertex::terminal_vertex(right_multiplier));
	  tmp.images_.insert(symbol, image_vertex);
		return tmp;
	}

	//! Returns the automorphism mapping #symbol to the product #left_multiplier #symbol.
	/**
	 * @param  left_multiplier  left multiplier
	 * @param symbol mapped terminal symbol, must be > 0
	 */
	static EndomorphismSLP left_multiplier(const TerminalSymbol& left_multiplier, const TerminalSymbol& symbol) {
	  assert(symbol > 0);
	  EndomorphismSLP tmp;
	  auto image_vertex = SLPVertex::concatenate(
      SLPVertex::terminal_vertex(left_multiplier),
      SLPVertex::terminal_vertex(symbol));
	  tmp.images_.insert(symbol, image_vertex);
		return tmp;
	}

	//! Returns the composition of endomorphisms specified by the range
	template<typename Iterator>
	static EndomorphismSLP composition(Iterator begin, Iterator end) {
	  return identity().compose_with(begin, end);
	}

	//! Returns the composition of #num endomorphisms produced by #producer
	/**
	 * @param num      the number of endomorphisms to compose
	 * @param producer endomorphisms generator (using operator())
	 */
  template<typename Generator>
  static EndomorphismSLP composition(unsigned int num, Generator& generator) {
    return identity().compose_with(num, generator);
  }

	//! Compose with the given endomorphism.
	EndomorphismSLP& operator*=(const EndomorphismSLP& a);

	//! Compose with the given endomorphism.
	EndomorphismSLP& operator*(const EndomorphismSLP& a) const {
		EndomorphismSLP result(*this);
		return result *= a;
	}

	//! Compose with endomorphisms specified by the range.
	template<typename Iterator>
  EndomorphismSLP& compose_with(Iterator begin, Iterator end) {
  for(;begin != end; ++begin)
    (*this) *= *begin;
  return *this;
	}

	//! Compose with #num automorphism produced by #producer
	/**
	   * @param num the number of endomorphisms to compose with
	   * @param producer endomorphisms generator (using operator())
	   */
  template<typename Generator>
  EndomorphismSLP& compose_with(unsigned int num, Generator& generator) {
    for (unsigned int i = 0; i < num; ++i)
      (*this) *= generator();
    return *this;
  }


	//! Returns the image of the terminal.
	SLPProducedWord image(const TerminalSymbol& t) const {
		return SLPProducedWord(slp(t));
	}

	//! Returns the root of the straight-line program representing the terminal symbol.
	SLPVertex slp(const TerminalSymbol& t) const {
	  const bool t_sign = t > 0;
	  const TerminalSymbol& positive_t = t_sign ? t : -t;
		auto result = images_.find(positive_t);
		if (result == images_.end()) //if it is not in images_, then it is the identity map.
			return SLPVertex::terminal_vertex(t);
		else
			return t_sign
			    ? result->second
			    : result->second.negate();
	}

private:
	//! The default constructor.
	EndomorphismSLP() {}

	EndomorphismSLP(const TerminalSymbol& inverted) {
		images_.insert(inverted,
						SLPVertex::terminal_vertex(inverted).negate());
	}

	//! The representation of images of terminal symbols by straight-line programs
	/**
	 * If there is no entry for a given terminal symbol, then its image is the terminal itself.
	 */
	std::map<TerminalSymbol, SLPVertex> images_;
};



//! Automorphisms generator
/**
 * @tparam TerminalSymbol          terminal symbol representation. We suppose it has a constructor taking index
 * @tparam RandomEngine            engine generating uniformly random non-negative numbers. See std::random library documentation.
 * @tparam TerminalSymbolIndexType terminal symbol index type
 */
template <typename TerminalSymbol = int,
  typename RandomEngine = std::default_random_engine,
  typename TerminalSymbolIndexType = int>
class UniformAutomorphismSLPGenerator {
public:
	typedef TerminalSymbolIndexType index_type;

	//! Constructs a generator of automorphisms of the free group of the given rank.
	/**
	 * @param rank free group rank > 0
	 */
	UniformAutomorphismSLPGenerator(index_type rank)
		: UniformAutomorphismSLPGenerator(rank, RandomEngine())
	{}

	//! Constructs a generator of automorphisms of the free group of the given rank.
	/**
	 * @param rank free group rank > 0
	 * @param seed random engine seed for creation of a new one
	 */
  UniformAutomorphismSLPGenerator(index_type rank, typename RandomEngine::result_type seed)
    : UniformAutomorphismSLPGenerator(rank, RandomEngine(seed))
  {}

	//! Constructs a generator of automorphisms of the free group of the given rank.
	/**
	 * @param rank free group rank > 0
	 * @param random_engine random engine
	 */
	UniformAutomorphismSLPGenerator(index_type rank, RandomEngine& random_engine)
    : RANK(rank),
      RIGHT_MULTIPLIERS_COUNT(2 * rank * (rank - 1)),
      INVERTERS_COUNT(rank),
      random_engine_(random_engine),
      random_distr_(0, COUNT - 1) {
	  assert(rank > 0);
  }

	//! Generates a random automorphism
	EndomorphismSLP<TerminalSymbol> operator()() {
	  index_type r_val = random_distr_(random_engine_);
	  if (r_val < MULTIPLIERS_COUNT) {
      const bool inverted = r_val % 2;
      r_val >>= 1;
      const index_type symbol_index = 1 + (r_val % RANK);
      const TerminalSymbol symbol(symbol_index);
      index_type multiplier_index = 1 + r_val / RANK;
      if (multiplier_index >= symbol_index)
        ++multiplier_index;//correction for the skipping of the symbol when pick the multiplier
      const TerminalSymbol multiplier(inverted ? - multiplier_index : multiplier_index);
      if (r_val < RIGHT_MULTIPLIERS_COUNT)
        return EndomorphismSLP<TerminalSymbol>::right_multiplier(symbol, multiplier);
      else
        return EndomorphismSLP<TerminalSymbol>::left_multiplier(multiplier, symbol);
    } else {
      return EndomorphismSLP<TerminalSymbol>::inverter(TerminalSymbol(r_val));
    }
	}


private:
	const index_type RANK;
	const index_type RIGHT_MULTIPLIERS_COUNT;
  const index_type LEFT_MULTIPLIERS_COUNT = RIGHT_MULTIPLIERS_COUNT;
  const index_type MULTIPLIERS_COUNT = RIGHT_MULTIPLIERS_COUNT + LEFT_MULTIPLIERS_COUNT;
  const index_type INVERTERS_COUNT;
  const index_type COUNT = MULTIPLIERS_COUNT + INVERTERS_COUNT;
  //! Random generator, which is a binding of provided distribution and uniform random engine.
  RandomEngine random_engine_;
  std::uniform_int_distribution<index_type> random_distr_;
};



template <typename TerminalSymbol>
EndomorphismSLP<TerminalSymbol>& EndomorphismSLP<TerminalSymbol>::operator*=(const EndomorphismSLP<TerminalSymbol>& a) {
	std::unordered_map<SLPVertex, SLPVertex> new_vertices;//a's vertices to new vertices correspondence

	for (auto iterator: a.images_) {
		const SLPVertex image_root = iterator->second;
		//for each root we go over the tree using postorder inspector,
		//attach terminal vertices to the roots of our endomorphism, and copy the tree above
		SLPPostorderInspector inspector(image_root);
		while (!inspector.inspection_ended()) {
			const SLPVertex& current_vertex = inspector.current_vertex();
			if (new_vertices.find(current_vertex) == new_vertices.end()) {//it was not copied yet
				if (current_vertex.is_terminal()) {//remap terminals to our roots
					auto our_root = slp(current_vertex.terminal_symbol());
					auto pair = std::make_pair(current_vertex,
							current_vertex.is_negative() ? our_root.negate() : our_root);
					new_vertices.insert(pair);
				} else {//for a nonterminal we already processed its children because postorder inspector
					auto left = new_vertices.find(current_vertex.left_child());
					auto right = new_vertices.find(current_vertex.right_child());
					auto new_vertex = SLPVertex::concatenate(left->second, right->second);
					new_vertices.insert(std::make_pair(current_vertex, new_vertex));
				}
			}
			inspector.go_to_next_vertex();
		}
	}

	//we update our roots to the new ones when necessary
	for (auto iterator: a.images_) {
		images_[iterator->first] = new_vertices.find(iterator->second);
	}
	return *this;
}



} /* namespace crag */
#endif /* CRAG_FREE_GROUPS_ENDOMORPHISM_SLP_H_ */
