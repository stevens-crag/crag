/*
 * FreeGroupAutomorhpism.cpp
 *
 *  Created on: Jan 29, 2013
 *      Author: pmorar
 */

#include "../include/FreeGroupAutomorhpism.h"

namespace crag {
FreeGroupAutomorphism::FreeGroupAutomorphism(size_type generators_num, TerminalSymbol inverted_terminal)
	: SLPSet(generators_num) {
	const int index = inverted_terminal - 1;
	roots[index] = roots[index].negate();
}

FreeGroupAutomorphism::FreeGroupAutomorphism(size_type generators_num,
		TerminalSymbol mapped_terminal,//TODO better name
		SLPVertex left_terminal_vertex,
		SLPVertex right_terminal_vertexl)
			: SLPSet(generators_num) {
	roots[mapped_terminal - 1] = SLPVertex::concatenate(left_terminal_vertex, right_terminal_vertexl);
}

FreeGroupAutomorphism& FreeGroupAutomorphism::composeWithNielsen(TerminalSymbol inverted_terminal) {
	const int index = inverted_terminal - 1;
	roots[index] = roots[index].negate();
	return this;
}

FreeGroupAutomorphism& FreeGroupAutomorphism::composeWithNielsen(TerminalSymbol mapped_terminal,
			SLPVertex left_terminal_vertex,
			SLPVertex right_terminal_vertex) {
	SLPVertex left = roots[left_terminal_vertex.terminal_symbol() - 1];
	if (left_terminal_vertex.is_negative())
		left = left.negate();

	SLPVertex right = roots[right_terminal_vertex.terminal_symbol() - 1];
		if (right_terminal_vertex.is_negative())
			right = right.negate();

	roots[mapped_terminal - 1] = SLPVertex::concatenate(left, right);

	return this;
}

FreeGroupAutomorphism& FreeGroupAutomorphism::apply(const FreeGroupAutomorphism& a) {
	if (generators_num() != a.generators_num())
		throw std::invalid_argument("Numbers of terminals are different");
	return *this;
}

} /* namespace crag */
