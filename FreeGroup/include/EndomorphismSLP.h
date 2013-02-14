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
template<typename TerminalSymbol = unsigned int>
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

	//! Compose automorphisms specified by the range
	template<typename Iterator>
	static EndomorphismSLP compose(Iterator begin, Iterator end) {
	  EndomorphismSLP tmp;
	  for(;begin != end; ++begin)
	    tmp *= *begin;
	  return tmp;
	}

	//! Compose #num automorphism produced by #producer
  template<typename Producer>
  static EndomorphismSLP compose(unsigned int num, Producer& producer) {
    EndomorphismSLP tmp;
    for (unsigned int i = 0; i < num; ++i)
      tmp *= producer();
    return tmp;
  }

	//! Compose with the given endomorphism.
	EndomorphismSLP& operator*=(const EndomorphismSLP& a);

	//! Compose with the given endomorphism.
	EndomorphismSLP& operator*(const EndomorphismSLP& a) const {
		EndomorphismSLP result(*this);
		return result *= a;
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
 * @tparam TerminalSymbol terminal symbol representation. We suppose that its default constructor creates 'null'
 * 	symbol, which we can increment using operator++ to enumerate first terminal symbols
 */
template <typename TerminalSymbol = unsigned int,
    typename RandomEngine = std::default_random_engine,
    typename RandomDistribution = std::uniform_int_distribution>
class UniformAutomorphismSLPGenerator {
public:
	typedef unsigned int size_type;//TODO it is TerminalSymbol index type
	typedef typename RandomEngine::result_type seed_type;

	//! Constructs a generator of automorphisms of the free group of the given rank.
	UniformAutomorphismSLPGenerator(size_type rank)
		: UniformAutomorphismSLPGenerator(rank, RandomEngine())
	{}

	//! Constructs a generator of automorphisms of the free group of the given rank.
	  UniformAutomorphismSLPGenerator(size_type rank, seed_type seed)
	    : UniformAutomorphismSLPGenerator(rank, RandomEngine(seed))
	  {}

	//! Constructs a generator of automorphisms of the free group of the given rank.
	UniformAutomorphismSLPGenerator(size_type rank, RandomEngine& random_engine)
    : RIGHT_MULTIPLIERS_COUNT(2 * rank * (rank - 1)),
      INVERTERS_COUNT(rank),
      random_(std::bind(RandomDistribution(0, COUNT), random_engine))
  {}

	//! Generates random automorphism as a composition of #num Nielsen automorphisms
	EndomorphismSLP<TerminalSymbol> operator()() {
	  return EndomorphismSLP<TerminalSymbol>::identity();//TODO make generation
	}


private:
	const size_type RIGHT_MULTIPLIERS_COUNT;
  const size_type LEFT_MULTIPLIERS_COUNT = RIGHT_MULTIPLIERS_COUNT;
  const size_type MULTIPLIERS_COUNT = RIGHT_MULTIPLIERS_COUNT + LEFT_MULTIPLIERS_COUNT;
  const size_type INVERTERS_COUNT;
  const size_type COUNT = MULTIPLIERS_COUNT + INVERTERS_COUNT;

  auto random_;
	//we don't use static fields because order of initialization is not defined
	//and static fields which are objects are not recommended to use
	std::vector<EndomorphismSLP<TerminalSymbol>> type1_nielsen;
	std::vector<EndomorphismSLP<TerminalSymbol>> type2_nielsen;
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
