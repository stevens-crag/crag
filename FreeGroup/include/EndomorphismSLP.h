/*
 * EndomorphsimSLP.h
 *
 *  Created on: Jan 29, 2013
 *      Author: pmorar
 */

#ifndef CRAG_FREE_GROUPS_INCLUDE_ENDOMORPHISM_SLP_H_
#define CRAG_FREE_GROUPS_INCLUDE_ENDOMORPHISM_SLP_H_

#include <map>
#include "SLPSet.h"

namespace crag {


/**
 * Represents a free group endomorphism using straight-line programs.
 * @tparam TerminalSymbol terminal symbols class
 */
template<typename TerminalSymbol>
class EndomorphismSLP {
public:

	//! Returns the identity automorphism
	static EndomorphismSLP identity() {
		return EndomorphismSLP();
	}

	//! Returns the Nielsen automorphism that inverts the specified terminal sybmbol.
	/**
	 * @param  inverted  inverted terminal symbol
	 */
	static EndomorphismSLP nielsenAutomorphism(TerminalSymbol inverted) {
		EndomorphismSLP tmp;
		tmp.images_.insert(inverted,
				SLPVertex::terminal_vertex(inverted).negate());
		return tmp;
	}

	//! Returns the Nielsen automorphism that maps a terminal to a product of two terminals
	/**
	 * The automorhism maps #multiplied terminal symbol to the product #multiplied #right-multiplier
	 * of two terminal symbols.
	 * @param	multiplied			mapped terminal symbol
	 * @param	right_multiplier	right multiplier
	 */
	static EndomorphismSLP nielsenAutomorphism(TerminalSymbol multiplied, TerminalSymbol right_multiplier) {
		auto image_vertex = SLPVertex::concatenate(
				SLPVertex::terminal_vertex(multiplied),
				SLPVertex::terminal_vertex(right_multiplier));
		EndomorphismSLP tmp;
		tmp.images_.insert(multiplied, image_vertex);
		return tmp;
	}

	//! Copy constructor
	explicit EndomorphismSLP(const EndomorphismSLP& e)
				: images_(e.images) {}

	//! Move constructor
	EndomorphismSLP(EndomorphismSLP&& e)
		: images_(std::move(e.images_)) {}

	//! Assignment operator
	EndomorphismSLP& operator=(const EndomorphismSLP& e) {
		if (this != &e) {
			images_ = e.images_;
		}
		return *this;
	}

	//! Move operator
	EndomorphismSLP& operator=(EndomorphismSLP&& e) {
		if (this != &e) {
			images_ = std::move(e.images_);
		}
		return *this;
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
		auto result = images_.find(t);
		if (result == images_.end()) //if it is not in images_, then it is the identity map.
			return SLPVertex::terminal_vertex(t);
		else
			return result->second;
	}

protected:
	//! The default constructor.
	EndomorphismSLP() {}
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
template <typename TerminalSymbol>
class AutomorphismSLPGenerator {
public:
	typedef typename std::vector<EndomorphismSLP<TerminalSymbol>>::size_type size_type;

	//! Interface of random integers generator.
	class RandomGeneratorInterface {//TODO(pmorar) check whether there are any interfaces like this in crag
		//! Returns a random integer from [0,.. n).
		size_type rnd(int n) = 0;
	};


	const size_type NIELSEN_TYPE_1_COUNT;
	const size_type NIELSEN_TYPE_2_COUNT;
	const size_type COUNT;

	//! Constructs a generator of automorphisms of the free group of the given rank.
	AutomorphismSLPGenerator(size_type rank)
		: NIELSEN_TYPE_1_COUNT(rank),
		  NIELSEN_TYPE_2_COUNT(rank * (rank - 1)),
		  COUNT(NIELSEN_TYPE_1_COUNT + NIELSEN_TYPE_2_COUNT),
		  type1_nielsen(NIELSEN_TYPE_1_COUNT),
		  type2_nielsen(NIELSEN_TYPE_2_COUNT)
	{
		TerminalSymbol terminal(++TerminalSymbol());//the first terminal symbol
		for (unsigned int i = 0; i < rank; ++i, ++terminal) {
			type1_nielsen.push_back(EndomorphismSLP<TerminalSymbol>::nielsenAutomorphism(terminal));
			TerminalSymbol right_multiplier_terminal(++TerminalSymbol());//the first terminal symbol
			for (unsigned int j = 0; j < rank; ++i, ++right_multiplier_terminal) {
				if (j == i)
					continue;
				type1_nielsen.push_back(terminal, right_multiplier_terminal);
			}
		}
	}

	//! Returns the Nielsen automorphism with the given index
	/**
	 * The indices from [0,..#NIELSEN_TYPE_1_COUNT) are for Nielsen automorphisms inverting a vertex,
	 * and the indices from [#NIELSEN_TYPE_1_COUNT,..#COUNT) are for Nielsen automorphism mapping
	 * a vertex to a product of another ones.
	 * @param index a number from [0,..,#COUNT)
	 */
	const EndomorphismSLP<TerminalSymbol>& get(size_type index) const {
			if (index >= NIELSEN_TYPE_1_COUNT) {
				return type2_nielsen[index - NIELSEN_TYPE_1_COUNT];
			}
			return type1_nielsen[index];
		}

	//! Returns the Nielsen automorphism inverting a vertex corresponding to the given index.
	/**
	 * @param index a number from [0,..,#NIELSEN_TYPE_1_COUNT)
	 */
	const EndomorphismSLP<TerminalSymbol>& getType1(size_type index) const {
		return type1_nielsen[index];
	}

	//! Returns the Nielsen automorphism mapping a vertex to the product of another ones corresponding to the given index.
	/**
	 * @param index a number from [0,..,#NIELSEN_TYPE_2_COUNT)
	 */
	const EndomorphismSLP<TerminalSymbol>& getType2(size_type index) const {
		return type2_nielsen[index];
	}

	//! Generates random automorphism as a composition of #num Nielsen automorphisms using random generator rnd
	EndomorphismSLP<TerminalSymbol> random_automorphism(unsigned int num, RandomGeneratorInterface* rnd) {
		auto endomorphism = EndomorphismSLP<TerminalSymbol>::identity();
		for (int i = 0; i < num; ++i) {
			endomorphism *= get(rnd->rnd(COUNT));
		}
		return endomorphism;
	}


private:
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
#endif /* CRAG_FREE_GROUPS_INCLUDE_ENDOMORPHISM_SLP_H_ */
