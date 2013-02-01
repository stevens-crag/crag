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

FreeGroupAutomorphism& FreeGroupAutomorphism::composeWith(const FreeGroupAutomorphism& a) {
	if (generators_num() != a.generators_num())
		throw std::invalid_argument("Numbers of terminals are different");
	return *this;
}

} /* namespace crag */
