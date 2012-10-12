/*
 * BaumslagGersten.h
 *
 *  Created on: 07.04.2011
 *      Author: Juern
 *
 * Solving the WP in the Baumslag-Gersten group
 * BG(1,2)=< a,t,b | a^t=a^2, a^b=t > (where
 * x^y = yxy^-1).
 *
 * There are two functions expecting different forms
 * of input. The first one expects the input to be
 * already in compressed form. The second one is more
 * classical and takes a mere string.
 */

#ifndef BAUMSLAGGERSTEN_H_
#define BAUMSLAGGERSTEN_H_

namespace BG
{

using namespace PC;

/**
 * This data type is used to represent the input of
 * the procedure solving the word problem in BG(1,2).
 * A word is given as a sequence of BGMonimials, i.e., letters b, b^-1 and
 * elements of the Baumslag-Solitar group
 * BS(1,2)=< a,t | a^t=a^2 >, where the latter are
 * represented by triple markings (U,X,K) in a common
 * power circuit, where:
 * i)	all U, X, and K must have disjoint supports,
 * ii)	all U must be sources,
 * iii)	all incoming arcs to X and K must originate
 * 		in the correspoding U, and
 * iv)	arcs from U to X must have the opposite sign
 * 		of the correspoding node-signs in X.
 */

struct BGMonomial
{
	enum {BSat, b}	type;
	//union
	//{
		Marking	U, X, K;	// if type == BSat
		int		expb;		// if type == b (must be +1 or -1)
	//};
};

bool solveWPinBG(PowerCircuit* pc, std::list<BGMonomial>& input);


/**
 * The input is a string consisting of the letters
 * a, b, t and their inverses A, B, T. An empty instance
 * of PowerCircuit is needed, too.
 */

bool solveWPinBG(PowerCircuit* pc, std::string input);

}

#endif /* BAUMSLAGGERSTEN_H_ */
