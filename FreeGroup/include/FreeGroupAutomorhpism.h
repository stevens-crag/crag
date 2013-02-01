/*
 * FreeGroupAutomorhpism.h
 *
 *  Created on: Jan 29, 2013
 *      Author: pmorar
 */

#ifndef FREEGROUPAUTOMORHPISM_H_
#define FREEGROUPAUTOMORHPISM_H_

#include "SLPSet.h"

namespace crag {



/**
 * Represents a free group automorphism.
 */
class FreeGroupAutomorphism: private crag::SLPSet {
public:
	typedef typename SLPSet::size_type size_type;

	FreeGroupAutomorphism() = delete;

	//! Constructs the identity automorphism of the free group with n generators
	FreeGroupAutomorphism(size_type generators_num): SLPSet(generators_num) {}

	//! Copy constructor
	explicit FreeGroupAutomorphism(const FreeGroupAutomorphism& x)
				: SLPSet(x) {}

	//! Move constructor
	FreeGroupAutomorphism(FreeGroupAutomorphism&& x)
		: SLPSet(std::move(x)) {}

	//! The Nielsen automorphism of the free group with n generators that just inverts the specified terminal and fixes all other terminals
	FreeGroupAutomorphism(size_type generators_num, TerminalSymbol inverted_terminal);

	//! The Nielsen automorphism of the free group with n generators that maps the specified terminal to the product of another two
	FreeGroupAutomorphism(size_type generators_num,
			TerminalSymbol mapped_terminal,//TODO better name
			SLPVertex left_terminal_vertex,
			SLPVertex right_terminal_vertexl);

	//! Composes the current automorphism with the specified Nielsen automorphism
	/**
	 * Changes the current automorphism so it represents the composition of itself with the given Nielsen
	 * automorphism, which inverts the specified terminal symbol.
	 *
	 * @param inverted_terminal the inverted terminal symbol
	 * @return itself after the composition
	 */
	FreeGroupAutomorphism& composeWithNielsen(TerminalSymbol inverted_terminal);

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
	FreeGroupAutomorphism& composeWithNielsen(TerminalSymbol mapped_terminal,
			SLPVertex left_terminal_vertex,
			SLPVertex right_terminal_vertexl);

	//! Applies the given automorphism to the current one and returns the resulting automorphism
	/**
	 * Applies the automorphism a to the given one. It constructs a new automorphism and returns it.
	 * The computational complexity is O(n ln n), where n is the size of the slp representing automorphism a.
	 */
	FreeGroupAutomorphism& apply(const FreeGroupAutomorphism& a);

	//! The image of the terminal t.
	/**
	 * @param t terminal (must be between 1 and generators_num())
	 */
	SLPProducedWord image(TerminalSymbol t) const {
		return SLPSet::produced_word(t - 1);
	}

	//! The root of SLP representing the image of the terminal t.
	/**
	 * @param t terminal (must be between 1 and generators_num())
	 */
	SLPVertex slp(TerminalSymbol t) const {
		return SLPSet::root(t - 1);
	}

	//! Returns the generators number of the free group
	size_type generators_num() const {
		return roots_num();
	}

};



} /* namespace crag */
#endif /* FREEGROUPAUTOMORHPISM_H_ */
