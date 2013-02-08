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
 * @tparam rank the rank of the free group
 */
template<unsigned int rank>
class EndomorphismSLP {
public:

	//! Returns the identity automorphism
	static EndomorphismSLP identity() {
		return EndomorphismSLP();
	}

	//! Copy constructor
	explicit EndomorphismSLP(const EndomorphismSLP<rank>& e)
				: images_(e.images) {}

	//! Move constructor
	EndomorphismSLP(EndomorphismSLP<rank>&& e)
		: images_(std::move(e.images_)) {}

	//! Assignment operator
	EndomorphismSLP<rank>& operator=(const EndomorphismSLP<rank>& e) {
		if (this != &e) {
			images_ = e.images_;
		}
		return *this;
	}

	//! Move operator
	EndomorphismSLP<rank>& operator=(EndomorphismSLP<rank>&& e) {
		if (this != &e) {
			images_ = std::move(e.images_);
		}
		return *this;
	}

	//! Compose with the given endomorphism.
	EndomorphismSLP<rank>& operator*=(const EndomorphismSLP<rank>& a);

	//! Compose with the given endomorphism.
	EndomorphismSLP<rank>& operator*(const EndomorphismSLP<rank>& a) const {
		EndomorphismSLP result(*this);
		return result *= a;
	}

	//! Returns the image of the terminal.
	SLPProducedWord image(const TerminalSymbol& t) const {
		return SLPProducedWord(slp(t));
	}

	//! Returns the root of the straight-line program representing the terminal sybmbol.
	SLPVertex slp(const TerminalSymbol& t) const {
		auto result = images_.find(t);
		if (result == images_.end())
			return SLPVertex::terminal_vertex(t);
		else
			return result->second;
	}

	//! Returns the generators number of the free group
	unsigned int get_rank() const {
		return rank;
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

template <unsigned int rank>
class NielsenAutomorphismSLP: public EndomorphismSLP<rank> {
public:
	//! Constructs the Nielsen automorphism that inverts the specified terminal sybmbol.
	/**
	 * @param  inverted  inverted terminal symbol
	 * @throws std::invalid_argument if inverted is out of bounds with respect to the rank
	 */
	NielsenAutomorphismSLP(TerminalSymbol inverted)
		throw(std::invalid_argument) {
			if (inverted > rank)
				throw std::invalid_argument();
			images_.insert(inverted,
					SLPVertex::terminal_vertex(inverted).negate());
	}

	//! Constructs the Nielsen automorphism that maps a terminal to a product of two terminals
	/**
	 * The automorhism maps #multiplied terminal symbol to the product #multiplied #right-multiplier
	 * of two terminal symbols.
	 * @param	multiplied			mapped terminal symbol
	 * @param	right_multiplier	right multiplier
	 * @throws std::invalid_argument if inverted or right_multiplier is out of bounds with respect to the rank
	 */
	NielsenAutomorphismSLP(TerminalSymbol multiplied,
		TerminalSymbol right_multiplier)
		throw(std::invalid_argument) {
			if (multiplied > rank || right_multiplier > rank)
					throw std::invalid_argument();
			auto image_vertex = SLPVertex::concatenate(
					SLPVertex::terminal_vertex(multiplied),
					SLPVertex::terminal_vertex(right_multiplier));
			images_.insert(multiplied, image_vertex);
		}

	bool is_nielsen() {
		if (images_.size() != 1 || images_.begin()->second.height() > 2)
			return false;
		return true;
	}
};


template <unsigned int rank>
class AutomorphismSLPGenerator {
public:
	typedef typename std::vector<NielsenAutomorphismSLP<rank>>::size_type size_type;

	//! Interface of random integers generator.
	class RandomGeneratorInterface {//TODO(pmorar) check whether there are any interfaces like this in crag
		//! Returns a random integer from [0,.. n).
		size_type rnd(int n) = 0;
	};


	static const unsigned int NIELSEN_TYPE_1_COUNT = rank;
	static const unsigned int NIELSEN_TYPE_2_COUNT = rank * (rank - 1);
	static const unsigned int COUNT = NIELSEN_TYPE_1_COUNT + NIELSEN_TYPE_2_COUNT;

	//! Constructor
	AutomorphismSLPGenerator()
		: type1_nielsen(NIELSEN_TYPE_1_COUNT),
		type2_nielsen(NIELSEN_TYPE_2_COUNT) {
		for (unsigned int i = 1; i <= rank; ++i) {
			type1_nielsen.push_back(NielsenAutomorphismSLP(TerminalSymbol(i)));
			for (unsigned int j = 1; j <= rank; ++i) {
				if (j == i)
					continue;
				type1_nielsen.push_back(NielsenAutomorphismSLP(TerminalSymbol(i), TerminalSymbol(j)));
			}
		}
	}

	const NielsenAutomorphismSLP<rank>& get(size_type index) const {
			if (index >= NIELSEN_TYPE_1_COUNT) {
				index -= NIELSEN_TYPE_1_COUNT;
				return type2_nielsen[index];
			}
			return type1_nielsen[index];
		}

	const NielsenAutomorphismSLP<rank>& getType1(size_type index) const {
		return type1_nielsen[index];
	}

	const NielsenAutomorphismSLP<rank>& getType2(size_type index) const {
		return type2_nielsen[index];
	}

	//! Generates random automorphism as a composition of #num Nielsen automorphisms using random generator rnd
	EndomorphismSLP<rank> random_automorphism(unsigned int num, RandomGeneratorInterface* rnd) {
		EndomorphismSLP<rank> endomorphism = EndomorphismSLP<rank>::identity();
		for (int i = 0; i < num; ++i) {
			endomorphism *= get(rnd->rnd(COUNT));
		}
		return endomorphism;
	}


private:
	//we don't use static fields because order of initialization is not defined
	//and static fields which are objects are not recommended to use
	std::vector<NielsenAutomorphismSLP<rank> > type1_nielsen;
	std::vector<NielsenAutomorphismSLP<rank> > type2_nielsen;
};



template <unsigned int rank>
EndomorphismSLP<rank>& EndomorphismSLP<rank>::operator*=(const EndomorphismSLP<rank>& a) {
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
