/*
 * EndomorphsimSLP.h
 *
 *  Created on: Jan 29, 2013
 *      Author: pmorar
 */

#ifndef ENDOMORPHISM_SLP_H_
#define ENDOMORPHISM_SLP_H_

#include "SLPSet.h"

namespace crag {



/**
 * Represents a free group endomorphism using straight-line programs.
 * @tparam rank the rank of the free group
 */
template<int rank>
class EndomorphismSLP {
public:
	typedef typename SLPSet::size_type size_type;

	//! Returns the identity automorphism
	const EndomorphismSLP<rank>& identity() {
		return _identity;
	}

	EndomorphismSLP& operator=(const EndomorphismSLP<rank>&) = delete;


	//! Copy constructor
	explicit EndomorphismSLP(const EndomorphismSLP<rank>& x)
				: slp_(x.slp_) {}

	//! Move constructor
	EndomorphismSLP(EndomorphismSLP<rank>&& x)
		: slp_(std::move(x.slp_)) {}

	//! The Nielsen automorphism that inverts the specified terminal and fixes all other terminals
	EndomorphismSLP(TerminalSymbol inverted_terminal)
	: slp_(generators_num) {
		slp_.invert_root(inverted_terminal - 1);
	}

	//! The Nielsen automorphism that maps the specified terminal to the product of another two
	EndomorphismSLP(TerminalSymbol mapped_terminal,
			SLPVertex left_terminal_vertex,
			SLPVertex right_terminal_vertex)
	: slp_(generators_num) {
		slp_.replace_root(mapped_terminal - 1,
						SLPVertex::concatenate(left_terminal_vertex, right_terminal_vertex));
	}

	//! Composes the current automorplhism with the specified Nielsen automorphism
	/**
	 * Changes the current automorphism so it represents the composition of itself with the given Nielsen
	 * automorphism, which inverts the specified terminal symbol.
	 *
	 * @param inverted_terminal the inverted terminal symbol
	 * @return itself after the composition
	 */
	EndomorphismSLP& composeWithNielsen(TerminalSymbol inverted_terminal) {
		slp_.invert_root(inverted_terminal - 1);
		return *this;
	}

	//! Composes the current automorphism with the specified Nielsen automorphism
	/**
	 * Changes the current automorphism so it represents the composition of itself with the given Nielsen
	 * automorphism, which maps the specified terminal to the product of another two.
	 *
	 * @param mapped_terminal terminal symbol being mapped
	 * @param left_terminal_vertex vertex representing left terminal in the product
	 * @param right_terminal_vertexlvertex representing right terminal in the product
	 * @return itself after the composition
	 */
	EndomorphismSLP& composeWithNielsen(TerminalSymbol mapped_terminal,
			SLPVertex left_terminal_vertex,
			SLPVertex right_terminal_vertex) {
		auto left = slp_.root(left_terminal_vertex.terminal_symbol() - 1);
		if (left_terminal_vertex.is_negative())
			left = left.negate();

		auto right = slp_.root(right_terminal_vertex.terminal_symbol() - 1);
		if (right_terminal_vertex.is_negative())
			right = right.negate();

		slp_.replace_root(mapped_terminal - 1,
				SLPVertex::concatenate(left, right));

		return *this;
	}

	//! Applies the given automorphism to the current one and returns the resulting automorphism
	/**
	 * Applies the automorphism a to the given one. It constructs a new automorphism and returns it.
	 * The computational complexity is O(n ln n), where n is the size of the slp representing automorphism a.
	 */
	EndomorphismSLP& apply(const EndomorphismSLP<rank>& a);

	//! The image of the terminal t.
	/**
	 * @param t terminal (must be between 1 and generators_num())
	 */
	SLPProducedWord image(TerminalSymbol t) const {
		return slp_.produced_word(t - 1);
	}

	//! The root of SLP representing the image of the terminal t.
	/**
	 * @param t terminal (must be between 1 and generators_num())
	 */
	SLPVertex slp(TerminalSymbol t) const {
		return slp_.root(t - 1);
	}

	//! Returns the generators number of the free group
	size_type generators_num() const {
		return rank;
	}

private:
	EndomorphismSLP(): slp_(rank) {}
	static const EndomorphismSLP<rank> _identity;
	SLPSet slp_;//representation

};


template <int rank>
EndomorphismSLP<rank>& EndomorphismSLP<rank>::apply(const EndomorphismSLP<rank>& a) {//TODO operator *
	std::unordered_map<SLPVertex, SLPVertex> vertices_map;

	{
		//We use postorder inspector to go over the tree and attach the terminals to
		//the roots of a and all the other vertices to them replicating the tree structure.
		//We do this differently for the first root because we have not copied any vertices yet.
		SLPPostorderInspector inspector(slp(1));
		while (!inspector.inspection_ended()) {
			const SLPVertex& current_vertex = inspector.current_vertex();
			if (current_vertex.is_terminal()) {//remap terminals to the roots of a
				auto a_vertex = a.slp(current_vertex.terminal_symbol());
				auto pair = std::make_pair(current_vertex,
						current_vertex.is_negative() ? a_vertex.negate() : a_vertex);
				vertices_map.insert(pair);
			} else {//for nonlterminal we already visited its children and copied them
				//because postorder
				auto left = vertices_map.find(current_vertex.left_child());
				auto right = vertices_map.find(current_vertex.right_child());
				auto u = SLPVertex::concatenate(left->second, right->second);
				auto pair = std::make_pair(current_vertex, u);
				vertices_map.insert(pair);
			}
			inspector.go_to_next_vertex();
		}
	}
	//now we have to move the other roots. For them part of the vertices is copied already
	//so we need to check for this condition
	for (TerminalSymbol terminal = 2, size = generators_num(); terminal <= size; ++terminal) {
		SLPPostorderInspector inspector(slp(terminal));//TODO make postorder iterator with function that tells when to stop
		while (!inspector.inspection_ended()) {
			const SLPVertex& current_vertex = inspector.current_vertex();
			if (vertices_map.find(current_vertex) == vertices_map.end()) {
				//we replace if it is not replaced already, this check may be removed
				//if the inspector does not go over such vertices
				if (current_vertex.is_terminal()) {//remap terminals to the roots of a
					auto a_vertex = a.slp(current_vertex.terminal_symbol());
					auto pair = std::make_pair(current_vertex,
							current_vertex.is_negative() ? a_vertex.negate() : a_vertex);
					vertices_map.insert(pair);
				} else {//for nonterminal we already visited its children and copied them
					//because postorder
					auto left = vertices_map.find(current_vertex.left_child());
					auto right = vertices_map.find(current_vertex.right_child());
					auto u = SLPVertex::concatenate(left->second, right->second);
					auto pair = std::make_pair(current_vertex, u);
					vertices_map.insert(pair);
				}
			}
			inspector.go_to_next_vertex();
		}
	}
	//replacing roots
	for (TerminalSymbol terminal = 1, n = generators_num(); terminal <= n; ++terminal) {
		auto v_iter = vertices_map.find(slp(terminal));
		slp_.replace_root(terminal - 1, v_iter->second);
	}
	return *this;
}




} /* namespace crag */
#endif /* ENDOMORPHISM_SLP_H_ */
