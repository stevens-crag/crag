/*
 * FreeGroupAutomorhpism.cpp
 *
 *  Created on: Jan 29, 2013
 *      Author: pmorar
 */

#include "../include/EndomorphismSLP.h"
#include <sstream>

namespace crag {
EndomorphismSLP::EndomorphismSLP(size_type generators_num, TerminalSymbol inverted_terminal)
	: slp_(generators_num) {
	slp_.invert_root(inverted_terminal - 1);
}

EndomorphismSLP::EndomorphismSLP(size_type generators_num,
		TerminalSymbol mapped_terminal,//TODO better name
		SLPVertex left_terminal_vertex,
		SLPVertex right_terminal_vertex)
			: slp_(generators_num) {
	slp_.replace_root(mapped_terminal - 1,
			SLPVertex::concatenate(left_terminal_vertex, right_terminal_vertex));
}

EndomorphismSLP& EndomorphismSLP::composeWithNielsen(TerminalSymbol inverted_terminal) {
	slp_.invert_root(inverted_terminal - 1);
	return *this;
}

EndomorphismSLP& EndomorphismSLP::composeWithNielsen(TerminalSymbol mapped_terminal,
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

EndomorphismSLP& EndomorphismSLP::apply(const EndomorphismSLP& a) {//TODO operator *
	if (generators_num() != a.generators_num()) {
		std::stringstream s("generators num == ");
		s << generators_num() << " != a.generators_num == " << a.generators_num();
		throw std::invalid_argument(s.str());
	}
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
