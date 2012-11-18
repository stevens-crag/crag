/*
 * PowerCircuit.cpp
 *
 *  Created on: 23.03.2011
 *      Author: Juern
 */

#include <assert.h>
#include <vector>
#include <math.h>
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <list>
#include <stdlib.h>
#include <stdarg.h>
#include <unistd.h>
#include "Sign.h"
#include "PowerCircuit.h"

namespace PC
{

//-----------------------------------------------
// Node


/**
 * Returns whether a node is reduced.
 */
bool Node::isReduced()
{
	return (pc->getRedNodeOrd(*this) != -1);
}

/**
 * Returns whether a node is defined.
 */
bool Node::isDefined()
{
	return (id != -1);
}

/**
 * Clones a node, that means a new node is inserted with a copy of the successor-marking.
 * The new node will be not reduced.
 */
Node Node::clone()
{
	return pc->cloneNode(*this);
}

/**
 *	Returns the index in the reduces part of the PC (thus the
 *	number of smaller nodes in the reduced part).
 *	If the node is not in the reduced part, -1 is returned.
 */
int Node::getOrder()
{
	int ord = pc->getRedNodeOrd(*this);
	assert(ord != -1);
	return ord;
}

/**
 * Returns a copy of the successor-marking.
 */
Marking Node::getSuccessorMarking()
{
	return pc->getSuccMarking(*this);
}

/**
 * Compares two node, if they are referencing to the same node of the same circuit.
 */
bool Node::operator==(Node n)
{
	return (id == n.id && pc == n.pc);
}

//-----------------------------------------------
// Marking


Marking::Marking(PowerCircuit* p, int i): pc(p), id(i)
{
	if(pc!=NULL&&id!=-1)
		pc->incMarkingRefCount(*this);
}

Marking::Marking(const Marking& m): pc(m.pc), id(m.id)
{
	if(pc!=NULL&&id!=-1)
		pc->incMarkingRefCount(*this);
}

Marking& Marking::operator= (const Marking& rhs)
{
	if(id != rhs.id || pc != rhs.pc)
	{
		if(pc!=NULL&&id!=-1)
			pc->decMarkingRefCount(*this);
		pc = rhs.pc;
		id = rhs.id;
		if(pc!=NULL&&id!=-1)
			pc->incMarkingRefCount(*this);
	}
	return *this;
}

Marking::~Marking()
{
	if(pc!=NULL&&id!=-1)
		pc->decMarkingRefCount(*this);
}

/**
 * Returns whether a marking is reduced.
 */
bool Marking::isReduced()
{
	return pc->isMarkingReduced(*this);
}

/**
 * Returns whether a marking is defined.
 */
bool Marking::isDefined()
{
	return (id != -1);
}

/**
 * Clones the marking. That means all nodes in the support of the marking are cloned and a new marking with
 * the same value as the original one (using only the newly created nodes) is created.
 */
Marking Marking::clone()
{
	return pc->cloneMarking(*this);
}

/**
 * Return a copy of the marking (that means a new marking with the same nodes).
 */
Marking Marking::copy()
{
	return pc->copyMarking(*this);
}

/**
 * Returns a list of all nodes in the support of the marking.
 */
std::list<Node> Marking::getNodes()
{
	return pc->getMarkingNodes(*this);
}

/**
 * Returns the sign of a given node n. 1 stands for PLUS, -1 for MINUS.
 */
int Marking::getSign(Node n)
{
	return pc->getNodeSignInMarking(n, *this);
}

/**
 * Returns a Marking increased by one. Does not change the marking.
 * The result is NOT reduced. If the marking is reduced, the result has support in the reduced
 * part of the circuit.
 */
Marking Marking::operator++(int dummy)
{
	return pc->incMarking(*this);
}

/**
 * Creates a new Marking with the inverse value of mark.
 */
Marking Marking::operator-()
{
	return pc->invMarking(*this);
}

/**
 * Creates a new Marking with the sum of the two markings. Make sure that for each node,
 * the sum lies in the range [-1,+1]!!!
 */
Marking Marking::operator+(const Marking& m)
{
	return pc->addMarkings(*this, m);
}

/**
 * Creates a new Marking r with the intersection of this and m, that is
 * r(p) = m(p), 	if (this(p)==m(p))
 * r(p) = 0 		else
 * for every node p.
 */
Marking Marking::operator&(const Marking& m)
{
	return pc->intersectMarkings(*this, m);
}

/**
 * Compares reduced Markings by their value.
 */
bool Marking::operator<(const Marking& m)
{
	return (pc->compareMarkings(*this,m)==-1);
}

/**
 * Compares reduced Markings by their value.
 */
bool Marking::operator>(const Marking& m)
{
	return (pc->compareMarkings(*this,m)==1);
}

/**
 * Compares reduced Markings by their value.
 */
bool Marking::operator<=(const Marking& m)
{
	return (pc->compareMarkings(*this,m)<=0);
}

/**
 * Compares reduced Markings by their value.
 */
bool Marking::operator>=(const Marking& m)
{
	return (pc->compareMarkings(*this,m)>=0);
}

/**
 * Compares reduced Markings by their value.
 */
bool Marking::operator==(const Marking& m)
{
	return (pc->compareMarkings(*this,m)==0);
}

/**
 * Compares reduced Markings by their value.
 */
bool Marking::operator!=(const Marking& m)
{
	return (pc->compareMarkings(*this,m)!=0);
}

/**
 * Returns whether the Marking m is a successor-marking.
 */
bool Marking::isSuccessorMarking()
{
	return pc->isSuccessorMarking(*this);
}

/**
 * If the Marking is a successor-marking,
 * the node whose successor-marking it is, is returned.
 * Otherwise the undefined node is returned.
 */
Node Marking::getIncidentNode()
{
	return pc->getIncidentNode(*this);
}

/**
 * Returns the smallest non-zero node of a reduced marking m.
 * If m is the zero-marking the undefined node is returned.
 * If m is not reduced the it might return another node than the smallest.
 */
Node Marking::getSmallestNode()
{
	assert(isReduced());
	return pc->getSmallestNodeInMarking(*this);
}



//-----------------------------------------------
// Power Circuit

Marking PowerCircuit::createMarkingFromNodes(unsigned int numNodes, ...)
{
	va_list val;
	va_start(val, numNodes);
	std::list<Node> nodeList;
	for(unsigned int i = 0; i < numNodes; i++)
	{
		Node n = va_arg(val, Node);
		nodeList.push_back(n);
	}
	va_end(val);
	return createMarking(nodeList);
}

void PowerCircuit::draw(std::string filename, Marking highlight1, Marking highlight2, Marking highlight3, Marking highlight4, Marking highlight5, Marking highlight6, Marking highlight7, Marking highlight8, Marking highlight9)
{
	std::ofstream of((filename + ".dot").c_str(), std::ios_base::out | std::ios_base::trunc);
	of << "digraph PC" << std::endl << "{" << std::endl;
	of << "graph [center=1];" << std::endl << std::endl;

	std::list<std::string> tmpfiles;
	std::list<Node> nodes = getNodes();
	for(std::list<Node>::iterator i = nodes.begin(); i != nodes.end(); i++)
	{
		unsigned int numMarkings = (highlight1.isDefined() && highlight1.getSign(*i) != ZERO ? 1 : 0) + (highlight2.isDefined() && highlight2.getSign(*i) != ZERO ? 1 : 0) + (highlight3.isDefined() && highlight3.getSign(*i) != ZERO ? 1 : 0) + (highlight4.isDefined() && highlight4.getSign(*i) != ZERO ? 1 : 0) + (highlight5.isDefined() && highlight5.getSign(*i) != ZERO ? 1 : 0);
		// Draw node
		std::stringstream ss; ss << i->id;
		of << i->id << " [shape=circle" + (numMarkings > 0 ? " image=\"" + filename + ss.str() + ".eps\"" : "") << " label=\"" << (i->isReduced() ? "R" :"") << "\"];" << std::endl;
		// Create the node's image according to markings
		Sign s;
		if(numMarkings > 0)
		{
			tmpfiles.push_back(filename + ss.str() + ".eps");
			std::ofstream of2((filename + ss.str() + ".eps").c_str(), std::ios_base::out | std::ios_base::trunc);
			of2 << "%!PS-Adobe-2.0 EPSF-2.0" << std::endl << "%%BoundingBox: 0 0 50 50" << std::endl << "%%EndComments"/* << std::endl << "/Times-Roman findfont 8 scalefont setfont"*/ << std::endl;
			int m = 0;
			if(highlight1.isDefined() && (s = highlight1.getSign(*i)) != ZERO)
			{
				of2 << "25 25 moveto" << std::endl << "25 25 25 " << (m * 360 / numMarkings) << " " << ((m + 1) * 360 / numMarkings) << " arc" << std::endl << "1 0 0 setrgbcolor fill" << std::endl;
				of2 << "14 /Times-Roman set_font" << std::endl << (int)(25.0f + 12.5f * cos(((float)m + 0.5) * 2.0f * 3.14159 / (float)numMarkings)) << " " << (int)(25.0f + 12.5f * sin(((float)m + 0.5) * 2.0f * 3.14159 / (float)numMarkings)) << " moveto" << std::endl << "0 0 0 setrgbcolor" << std::endl << (s == PLUS ? "(+) show" : "(-) show") << std::endl;
				m++;
			}
			if(highlight2.isDefined() && (s = highlight2.getSign(*i)) != ZERO)
			{
				of2 << "25 25 moveto" << std::endl << "25 25 25 " << (m * 360 / numMarkings) << " " << (int)((m + 1) * 360 / numMarkings) << " arc" << std::endl << "0.7 0 0 setrgbcolor fill" << std::endl;
				of2 << "14 /Times-Roman set_font" << std::endl << (int)(25.0f + 12.5f * cos(((float)m + 0.5) * 2.0f * 3.14159 / (float)numMarkings)) << " " << (int)(25.0f + 12.5f * sin(((float)m + 0.5) * 2.0f * 3.14159 / (float)numMarkings)) << " moveto" << std::endl << "0 0 0 setcolor" << std::endl << (s == PLUS ? "(+) show" : "(-) show") << std::endl;
				m++;
			}
			if(highlight3.isDefined() && (s = highlight3.getSign(*i)) != ZERO)
			{
				of2 << "25 25 moveto" << std::endl << "25 25 25 " << (m * 360 / numMarkings) << " " << (int)((m + 1) * 360 / numMarkings) << " arc" << std::endl << "0.4 0 0 setrgbcolor fill" << std::endl;
				of2 << "14 /Times-Roman set_font" << std::endl << (int)(25.0f + 12.5f * cos(((float)m + 0.5) * 2.0f * 3.14159 / (float)numMarkings)) << " " << (int)(25.0f + 12.5f * sin(((float)m + 0.5) * 2.0f * 3.14159 / (float)numMarkings)) << " moveto" << std::endl << "0 0 0 setcolor" << std::endl << (s == PLUS ? "(+) show" : "(-) show") << std::endl;
				m++;
			}
			if(highlight4.isDefined() && (s = highlight4.getSign(*i)) != ZERO)
			{
				of2 << "25 25 moveto" << std::endl << "25 25 25 " << (m * 360 / numMarkings) << " " << (int)((m + 1) * 360 / numMarkings) << " arc" << std::endl << "0 1 0 setrgbcolor fill" << std::endl;
				of2 << "14 /Times-Roman set_font" << std::endl << (int)(25.0f + 12.5f * cos(((float)m + 0.5) * 2.0f * 3.14159 / (float)numMarkings)) << " " << (int)(25.0f + 12.5f * sin(((float)m + 0.5) * 2.0f * 3.14159 / (float)numMarkings)) << " moveto" << std::endl << "0 0 0 setcolor" << std::endl << (s == PLUS ? "(+) show" : "(-) show") << std::endl;
				m++;
			}
			if(highlight5.isDefined() && (s = highlight5.getSign(*i)) != ZERO)
			{
				of2 << "25 25 moveto" << std::endl << "25 25 25 " << (m * 360 / numMarkings) << " " << (int)((m + 1) * 360 / numMarkings) << " arc" << std::endl << "0 0.7 0 setrgbcolor fill" << std::endl;
				of2 << "14 /Times-Roman set_font" << std::endl << (int)(25.0f + 12.5f * cos(((float)m + 0.5) * 2.0f * 3.14159 / (float)numMarkings)) << " " << (int)(25.0f + 12.5f * sin(((float)m + 0.5) * 2.0f * 3.14159 / (float)numMarkings)) << " moveto" << std::endl << "0 0 0 setcolor" << std::endl << (s == PLUS ? "(+) show" : "(-) show") << std::endl;
				m++;
			}
			if(highlight6.isDefined() && (s = highlight6.getSign(*i)) != ZERO)
			{
				of2 << "25 25 moveto" << std::endl << "25 25 25 " << (m * 360 / numMarkings) << " " << (int)((m + 1) * 360 / numMarkings) << " arc" << std::endl << "0 0.4 0 setrgbcolor fill" << std::endl;
				of2 << "14 /Times-Roman set_font" << std::endl << (int)(25.0f + 12.5f * cos(((float)m + 0.5) * 2.0f * 3.14159 / (float)numMarkings)) << " " << (int)(25.0f + 12.5f * sin(((float)m + 0.5) * 2.0f * 3.14159 / (float)numMarkings)) << " moveto" << std::endl << "0 0 0 setcolor" << std::endl << (s == PLUS ? "(+) show" : "(-) show") << std::endl;
				m++;
			}
			if(highlight7.isDefined() && (s = highlight7.getSign(*i)) != ZERO)
			{
				of2 << "25 25 moveto" << std::endl << "25 25 25 " << (m * 360 / numMarkings) << " " << (int)((m + 1) * 360 / numMarkings) << " arc" << std::endl << "0 0 1 setrgbcolor fill" << std::endl;
				of2 << "14 /Times-Roman set_font" << std::endl << (int)(25.0f + 12.5f * cos(((float)m + 0.5) * 2.0f * 3.14159 / (float)numMarkings)) << " " << (int)(25.0f + 12.5f * sin(((float)m + 0.5) * 2.0f * 3.14159 / (float)numMarkings)) << " moveto" << std::endl << "0 0 0 setcolor" << std::endl << (s == PLUS ? "(+) show" : "(-) show") << std::endl;
				m++;
			}
			if(highlight8.isDefined() && (s = highlight8.getSign(*i)) != ZERO)
			{
				of2 << "25 25 moveto" << std::endl << "25 25 25 " << (m * 360 / numMarkings) << " " << (int)((m + 1) * 360 / numMarkings) << " arc" << std::endl << "0 0 0.7 setrgbcolor fill" << std::endl;
				of2 << "14 /Times-Roman set_font" << std::endl << (int)(25.0f + 12.5f * cos(((float)m + 0.5) * 2.0f * 3.14159 / (float)numMarkings)) << " " << (int)(25.0f + 12.5f * sin(((float)m + 0.5) * 2.0f * 3.14159 / (float)numMarkings)) << " moveto" << std::endl << "0 0 0 setcolor" << std::endl << (s == PLUS ? "(+) show" : "(-) show") << std::endl;
				m++;
			}
			if(highlight9.isDefined() && (s = highlight9.getSign(*i)) != ZERO)
			{
				of2 << "25 25 moveto" << std::endl << "25 25 25 " << (m * 360 / numMarkings) << " " << (int)((m + 1) * 360 / numMarkings) << " arc" << std::endl << "0 0 0.4 setrgbcolor fill" << std::endl;
				of2 << "14 /Times-Roman set_font" << std::endl << (int)(25.0f + 12.5f * cos(((float)m + 0.5) * 2.0f * 3.14159 / (float)numMarkings)) << " " << (int)(25.0f + 12.5f * sin(((float)m + 0.5) * 2.0f * 3.14159 / (float)numMarkings)) << " moveto" << std::endl << "0 0 0 setcolor" << std::endl << (s == PLUS ? "(+) show" : "(-) show") << std::endl;
				m++;
			}
			of2 << "%%EOF";
			of2.close();
		}
		// Draw arcs to other nodes
		Marking succMark = i->getSuccessorMarking();
		std::list<Node> succNodes = succMark.getNodes();
		for(std::list<Node>::iterator j = succNodes.begin(); j != succNodes.end(); j++)
		{
			s = succMark.getSign(*j);
			if(s != ZERO)
				of << i->id << " -> " << j->id << " [label=\"" << (s == PLUS ? "+" : "-") << "\"];" << std::endl;
		}
	}
	of << "}" << std::endl;
	of.close();
	// create image
	system(("dot -Tps -o" + filename + ".ps " + filename + ".dot").c_str());
	// remove temporary files
	unlink((filename + ".dot").c_str());
	while(!tmpfiles.empty())
	{
		unlink(tmpfiles.back().c_str());
		tmpfiles.pop_back();
	}
}

int PowerCircuit::getNumEdges()
{
	int result = 0;
	std::list<Node> nodeList = getNodes();
	for(std::list<Node>::iterator nodeIter= nodeList.begin(); nodeIter!= nodeList.end(); nodeIter++)
		result += nodeIter->getSuccessorMarking().getNodes().size();
	return result;
}

int PowerCircuit::getNumNodes()
{
	return getNodes().size();
}

int PowerCircuit::getNumMarkings()
{
	return getMarkings().size();
}

void PowerCircuit::printStatistics(std::ostream& os)
{
	os << "In this PC there are " << getNumMarkings() << " markings and " << getNumNodes() << " nodes." << std::endl;
	os << "Edges in underlying graph: " << getNumEdges() << ". " << std::endl << std::endl;
}

}
