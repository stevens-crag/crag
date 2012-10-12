/*
 * BaumslagGersten.cpp
 *
 *  Created on: 07.04.2011
 *      Author: Juern
 */

#include <string>
#include <cstring>
#include <vector>
#include <iostream>
#include <iomanip>
#include <list>
#include <assert.h>
#include <math.h>
#include "Sign.h"
#include "PowerCircuit.h"
#include "BaumslagGersten.h"
#include  "SignMatrix.h"
#include "PowerCircuitCompMatrix.h"

namespace BG
{

/**
 * Multiplies two triple markings in a sequence (the ones
 * at pos1 and pos2). The result is written to pos1 and
 * pos2 is removed, thus decreasing the length by one.
 * It is NOT checked, whether the two elements actually
 * exist and are both of type BSat!
 */
void multiply(PowerCircuit* pc, std::list<BGMonomial>& seq, std::list<BGMonomial>::iterator pos1, std::list<BGMonomial>::iterator& pos2)
{
	pc->connectInv(pos1->U, pos2->X);
	pc->connect(pos2->U, pos1->K);
	pos1->U = pos1->U + pos2->U;
	pos1->X = pos1->X + pos2->X;
	pos1->K = pos1->K + pos2->K;
	pos2 = seq.erase(pos2);
}

/**
 * See comment in BaumslagGersten.h
 */
bool solveWPinBG(PowerCircuit* pc, std::list<BGMonomial>& input)
{
	// Exclude trivial case
	if(input.size() == 0)
		return true;

	// Find monomials of type BS(a,t) standing next to each other
	// and multiply them.
	if(input.size() > 1)
	{
		std::list<BGMonomial>::iterator current = input.begin(), next = input.begin(); next++;
		while(next != input.end())
		{
			if(current->type == BGMonomial::BSat && next->type == BGMonomial::BSat)
				multiply(pc, input, current, next);
			else
			{
				current++;
				next++;
			}
		}
	}

	// Perform left-to-right Britton reductions:
	// Look for bb^-1 or b^-1b or b(u,x,k)b^-1 or b^-1(u,x,k)b
	// and (if possible) reduce.
	if(input.size() > 1)
	{
		std::list<BGMonomial>::iterator current = input.begin(), next = input.begin(); next++;
		while(next != input.end())
		{
			if(current->type == BGMonomial::b)
			{
				std::list<BGMonomial>::iterator secNext = next; secNext++;
				if(next->type == BGMonomial::b) // next letter is a b or b^-1
				{
					if(current->expb * next->expb == +1) // bb or b^-1b^-1
					{
						current++;
						next++;
					}
					else // bb^-1 or b^-1b
						next = input.erase(current, next);
				}
				else if(secNext != input.end() && secNext->expb * current->expb == -1) // b(U,X,K)b^-1 or b^-1(U,X,K)b
				{
					// Reduce a copy of the whole power circuit
					// while retaining the markings U, -X, and X+K.
					std::vector<Marking> markings; markings.resize(4);
					markings[0] = pc->createMarking(0); markings[1] = next->U; markings[2] = -next->X; markings[3] = next->X + next->K;
					PowerCircuit* pcClone = pc->clone(markings);
					pcClone->reduce();

					bool brittonReductionOccured = false;
					if(current->expb == +1) // b(U,X,K)b^-1
					{
						if(markings[3] == markings[0]) // X + K == 0
						{
							std::list<Node> succ; succ.push_back(markings[1].getSmallestNode());
							std::vector<Marking> markings2; markings2.resize(2);
							markings2[0] = pcClone->createMarking(succ); markings2[1] = markings[2];
							PowerCircuit* pcClone2 = pcClone->clone(markings2);
							pcClone2->reduce();
							if(markings2[0] >= markings2[1]) // U*2^X is integer
							{
								// multiply U by 2^X
								pc->connect(next->U, next->X);
								// set X:=U or K:= U depending on the sign of U
								if(markings[1] >= markings[0]) // U >= 0
								{
									next->X = pc->createMarking(0);
									next->K = next->U;
								}
								else
								{
									next->X = next->U;
									next->K = pc->createMarking(0);
								}
								next->U = pc->createMarking(0);
								brittonReductionOccured = true;
							}
//							delete pcClone2;
						}
					}
					else // b^-1(U,X,K)b
					{
						if(markings[1] == markings[0]) // U == 0
						{
							// remove U from the power circuit (thus making X and K sources)
							pc->remove(next->U);
							// make X + K the new U
							next->U = next->X + next->K;
							next->X = pc->createMarking(0);
							next->K = pc->createMarking(0);
							brittonReductionOccured = true;
						}
					}
//					delete pcClone;

					// So far, if there was a Britton reduction, we have only done
					// the swap. Now, remove the two surrounding b^(+/-1).
					if(brittonReductionOccured == true)
					{
						// Erase the b and b^-1 from the input sequence and,
						// if necessary, multiply triple markings.
						input.erase(secNext); input.erase(current);
						if(next == input.begin())
						{
							current = next; next++;
							if(next != input.end() && next->type == BGMonomial::BSat)
								multiply(pc, input, current, next);
						}
						else
						{
							current = next; current--;
							if(next != input.end())
							{
								secNext = next; secNext++;
								if(secNext != input.end() && secNext->type == BGMonomial::BSat)
									multiply(pc, input, next, secNext);
							}
							if(current->type == BGMonomial::BSat)
							{
								multiply(pc, input, current, next);
								if(current != input.begin())
								{
									current--;
									next--;
								}
							}
						}
					}
				}
				else
				{
					current++;
					next++;
				}
			}
			else
			{
				current++;
				next++;
			}
		}
	}
	// We now have a Britton-reduced sequence. However, it might
	// still contain elements of type BS(1,2) that are equal to 1.
	// In a final step remove these:
	std::list<Marking> markings;
	std::list<BGMonomial>::iterator iter = input.begin();
	while(iter != input.end())
	{
		if(iter->type == BGMonomial::BSat)
		{
			markings.push_back(iter->U);
			markings.push_back(iter->X + iter->K);
			iter++;
		}
	}
	Marking zero = pc->createMarking(0);
	pc->reduce();
	std::list<Marking>::iterator mIter = markings.begin();
	iter = input.begin();
	while(iter != input.end())
	{
		if(iter->type == BGMonomial::BSat)
		{
			if(*(mIter++) == zero && *(mIter++) == zero)
				iter = input.erase(iter);
			else
				iter++;
		}
	}
	return (input.size() == 0);
}

/*
 * The following two procedures deviate from the theoretically
 * best algorithm by not reducing just copies of the circuit but
 * the circuit itself. As the invariants about incoming and
 * outgoing edges of the markings U, X, and K cannot be kept up,
 * cloning becomes necessary before any operation + and connect.
 */

void multiplyWithCloning(PowerCircuit* pc, std::list<BGMonomial>& seq, std::list<BGMonomial>::iterator pos1, std::list<BGMonomial>::iterator& pos2)
{
	pos1->U = pos1->U.clone();
	pos2->X = pos2->X.clone();
	pc->connectInv(pos1->U, pos2->X);
	pos2->U = pos2->U.clone();
	pos1->K = pos1->K.clone();
	pc->connect(pos2->U, pos1->K);
	pos1->U = pos1->U + pos2->U;
	pos1->X = pos1->X + pos2->X.clone();
	pos1->K = pos1->K + pos2->K;
	pos2 = seq.erase(pos2);
}

bool solveWPinBGKeepingReductions(PowerCircuit* pc, std::list<BGMonomial>& input)
{
	// Exclude trivial case
	if(input.size() == 0)
		return true;

	pc->draw("begin");
	// Find monomials of type BS(a,t) standing next to each other
	// and multiply them.
	if(input.size() > 1)
	{
		std::list<BGMonomial>::iterator current = input.begin(), next = input.begin(); next++;
		while(next != input.end())
		{
			if(current->type == BGMonomial::BSat && next->type == BGMonomial::BSat)
				multiplyWithCloning(pc, input, current, next);
			else
			{
				current++;
				next++;
			}
		}
	}

	// Perform left-to-right Britton reductions:
	// Look for bb^-1 or b^-1b or b(u,x,k)b^-1 or b^-1(u,x,k)b
	// and (if possible) reduce.
	if(input.size() > 1)
	{
		std::list<BGMonomial>::iterator current = input.begin(), next = input.begin(); next++;
		while(next != input.end())
		{
			if(current->type == BGMonomial::b)
			{
				std::list<BGMonomial>::iterator secNext = next; secNext++;
				if(next->type == BGMonomial::b) // next letter is a b or b^-1
				{
					if(current->expb * next->expb == +1) // bb or b^-1b^-1
					{
						current++;
						next++;
					}
					else // bb^-1 or b^-1b
					{
						next = input.erase(current, next);
						current = next; current--;
					}
				}
				else if(secNext != input.end() && secNext->expb * current->expb == -1) // b(U,X,K)b^-1 or b^-1(U,X,K)b
				{

					// Reduce a copy of the whole power circuit
					// while retaining the markings U, -X, and X+K.
					std::vector<Marking> markings; markings.resize(4);
					markings[0] = pc->createMarking(0); markings[1] = next->U; markings[2] = -next->X; markings[3] = next->X + next->K;
					PowerCircuit* pcClone = pc;
					pcClone->reduce();


					bool brittonReductionOccured = false;
					if(current->expb == +1) // b(U,X,K)b^-1
					{

						if(markings[3] == markings[0]) // X + K == 0
						{


							std::list<Node> succ; succ.push_back(markings[1].getSmallestNode());
							std::vector<Marking> markings2; markings2.resize(2);
							markings2[0] = pcClone->createMarking(succ); markings2[1] = markings[2];
							PowerCircuit* pcClone2 = pc;
							pcClone2->reduce();
							if(markings2[0] >= markings2[1]) // U*2^X is integer
							{


								// multiply U by 2^X
								next->U = next->U.clone();
								next->X = next->X.clone();
								pc->connect(next->U, next->X);

								// set X:=U or K:= U depending on the sign of U
								if(markings[1] >= markings[0]) // U >= 0
								{
									next->X = pc->createMarking(0);
									next->K = next->U;
								}
								else
								{
									next->X = next->U;
									next->K = pc->createMarking(0);
								}
								next->U = pc->createMarking(0);
								brittonReductionOccured = true;

							}
							else
							{
								current++;
								next++;
							}

						}
						else
						{
							current++;
							next++;
						}
					}
					else // b^-1(U,X,K)b
					{
						if(markings[1] == markings[0]) // U == 0
						{
							// remove U from the power circuit (thus making X and K sources)
							pc->remove(next->U);
							// make X + K the new U
							next->U = next->X + next->K;
							next->X = pc->createMarking(0);
							next->K = pc->createMarking(0);
							brittonReductionOccured = true;
						}
						else
						{
							current++;
							next++;
						}
					}

					// So far, if there was a Britton reduction, we have only done
					// the swap. Now, remove the two surrounding b^(+/-1).
					if(brittonReductionOccured == true)
					{
						// Erase the b and b^-1 from the input sequence and,
						// if necessary, multiply triple markings.
						input.erase(secNext); input.erase(current);
						if(next == input.begin())
						{
							current = next; next++;
							if(next != input.end() && next->type == BGMonomial::BSat)
								multiplyWithCloning(pc, input, current, next);
						}
						else
						{
							current = next; current--;
							if(next != input.end())
							{
								secNext = next; secNext++;
								if(secNext != input.end() && secNext->type == BGMonomial::BSat)
									multiplyWithCloning(pc, input, next, secNext);
							}
							if(current->type == BGMonomial::BSat)
							{
								multiplyWithCloning(pc, input, current, next);
								if(current != input.begin())
								{
									current--;
									next--;
								}
							}
						}
					}
				}
				else
				{
					current++;
					next++;
				}
			}
			else
			{
				current++;
				next++;
			}
		}
	}

	// We now have a Britton-reduced sequence. However, it might
	// still contain elements of type BS(1,2) that are equal to 1.
	// In a final step remove these:
	std::list<Marking> markings;
	std::list<BGMonomial>::iterator iter = input.begin();
	while(iter != input.end())
	{
		if(iter->type == BGMonomial::BSat)
		{
			markings.push_back(iter->U);
			markings.push_back(iter->X + iter->K);
		}
		iter++;
	}
	Marking zero = pc->createMarking(0);
	pc->reduce();

	std::list<Marking>::iterator mIter = markings.begin();
	iter = input.begin();
	while(iter != input.end())
	{
		if(iter->type == BGMonomial::BSat)
		{
			if(*(mIter++) == zero && *(mIter++) == zero)
				iter = input.erase(iter);
			else
				iter++;
		}
		else
			iter++;
	}

	return (input.size() == 0);
}

/**
 * See comment in BaumslagGersten.h
 */
bool solveWPinBG(PowerCircuit* pc, std::string input)
{
	// Create a power circuit and a lot of disjoint +1 and -1 markings
	Marking						zero = pc->createMarking(0);
	Marking						one = pc->createMarking(1);
	BGMonomial					bgMon;
	std::list<BGMonomial>		sequence;
	for(unsigned int i = 0; i < input.size(); i++)
	{
		if(input[i] == ' ' || input[i] == '\t' || input[i] == '(' || input[i] == ')') // ignore whitespace and parantheses
			continue;
		if(input[i] == 'a' || input[i] == 'A')
		{
			bgMon.type = BGMonomial::BSat;
			bgMon.U = (input[i] == 'a' ? one.clone() : -one.clone());
			bgMon.X = zero;
			bgMon.K = zero;
		}
		else if(input[i] == 't' || input[i] == 'T')
		{
			bgMon.type = BGMonomial::BSat;
			bgMon.U = zero;
			bgMon.X = (input[i] == 't' ? zero : -one.clone());
			bgMon.K = (input[i] == 't' ? one.clone() : zero);
		}
		else if(input[i] == 'b' || input[i] == 'B')
		{
			bgMon.type = BGMonomial::b;
			bgMon.expb = (input[i] == 'b' ? +1 : -1);
		}
		else
			assert(false);
		sequence.push_back(bgMon);
	}


	//return solveWPinBG(pc, sequence);
	return solveWPinBGKeepingReductions(pc, sequence);
}

}
