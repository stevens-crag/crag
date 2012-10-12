/*
 * PowerCircuitGraph.cpp
 *
 *  Created on: 23.03.2011
 *      Author: ArminW
 */
#include <assert.h>
#include <vector>
#include <math.h>
#include <set>
#include <string>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <list>
#include <stdlib.h>
#include <map>
#include "Sign.h"
#include "PowerCircuit.h"
#include "PowerCircuitGraph.h"

namespace PC
{

PowerCircuitGraph::PowerCircuitGraph()
	:
	numNodes(0),
	numMarkings(0),
	numReducedNodes(0),
	firstDeletedNode(-1),
	firstDeletedMarking(-1)
{
}


PowerCircuitGraph::~PowerCircuitGraph()
{
}


//-------------------------------------------------------------------------------------------
// Functions for checking integrity of nodes and markings

/**
 * checks if a marking is valid
 */
bool PowerCircuitGraph::checkMarkingValid(const Marking& m) const
{
	if((unsigned int)m.id >= markings.size() || m.id < 0)
		return false;
	if(markings[m.id].deleted)
		return false;
	if(m.pc != this)
		return false;
	return true;
}


/**
 *  checks if a node is valid
 */
bool PowerCircuitGraph::checkNodeValid(Node n) const
{
	if((unsigned int)n.id >= nodes.size() || n.id < 0)
		return false;
	if(nodes[n.id].deleted)
		return false;
	if(n.pc != this)
		return false;
	return true;
}

//-------------------------------------------------------------------------------------------
// Functions to delete nodes and markings

/**
 * Deletes a node, its lambda-marking and also its corresponding row in the matrix.
 * If a reduced node is deleted, the user has to assure, that this does not change the property of being reduced
 * for the other reduced nodes.
 * Note that the values of other nodes and markings might be changed due to deletion of this node.
 */
void PowerCircuitGraph::deleteNode(Node node)
{
	if(nodes[node.id].indexInOrder < numReducedNodes)
	{
		std::cout<<"Attention: Delete Reduced Node!!"<<std::endl;
		if(nodes[node.id].indexInOrder  > 0)
			nodes[nodeOrder[nodes[node.id].indexInOrder - 1].id].BV = false;

		// no nodes are reduced anymore (except the one-node)
		for(unsigned int i = 0 ; i < numReducedNodes ; i++)
		{
			nodes[nodeOrder[i].id].BV = false;
		}
		numReducedNodes = (nodes[node.id].indexInOrder==0)?0:1;
	}

	for(unsigned int i = nodes[node.id].indexInOrder ; i < numNodes - 1 ; i++)
	{
		nodeOrder[i] = nodeOrder[i+1];
		nodes[nodeOrder[i].id].indexInOrder = i;
	}
	nodeOrder.pop_back();

	nodes[node.id].deleted = true;
	nodes[node.id].successors.clear();
	nodes[node.id].nextDeleted = firstDeletedNode;
	firstDeletedNode = node.id;
	numNodes--;
	assert(numNodes == nodeOrder.size());
}

/**
 * Deletes the marking mark.
 * Does NOT delete nodes which are not needed anymore.
 */
void PowerCircuitGraph::deleteMarking(Marking& mark)
{
	assert(checkMarkingValid(mark));

	markings[mark.id].nodes.clear();
	markings[mark.id].deleted = true;
	markings[mark.id].reduced = false;
	markings[mark.id].nextDeleted = firstDeletedMarking;
	numMarkings--;
	firstDeletedMarking = mark.id;
	mark.id = -1;
	mark.pc = NULL;
}


//-------------------------------------------------------------------------------------------
// Internal functions to create nodes and markings

/**
 * Creates a new marking out of nodeList.
 * Markings internally MUST NOT be created through other methods than this!!!
 */
Marking PowerCircuitGraph::newMarking(std::list< std::pair<Node,Sign> >  nodeList = std::list< std::pair<Node,Sign> >() )
{
	int n;
	if(firstDeletedMarking == -1)
	{
		markings.resize(numMarkings + 1);
		n = numMarkings;
	}
	else
	{
		n = firstDeletedMarking;
		firstDeletedMarking = markings[n].nextDeleted;
	}
	markings[n].deleted = false;
	markings[n].reduced = false;
	markings[n].refCount = 0;
	markings[n].nodes = nodeList;

	numMarkings++;
	return Marking(this, n);
}


/**
 * Takes list of signed nodes and turns it into the successor-marking of a new node.
 * The newly created node is returned.
 * Nodes internally SHOULD NOT be created through other methods than this (except for cloning PCs)!!!
 */
Node PowerCircuitGraph::newNode(std::list< std::pair<Node,Sign> >  nodeList = std::list< std::pair<Node,Sign> >() )
{
	int nodeNum;

	// insert node into static list of all nodes
	if(firstDeletedNode == -1)
	{
		nodes.resize(numNodes + 1);
		nodeNum = numNodes;
	}
	else
	{
		nodeNum = firstDeletedNode;
		firstDeletedNode = nodes[nodeNum].nextDeleted;
	}

	nodes[nodeNum].deleted = false;
	nodes[nodeNum].successors = nodeList;

	nodes[nodeNum].BV = false;
	nodes[nodeNum].indexInOrder = numNodes;

	assert(nodeOrder.size() == numNodes);
	nodeOrder.push_back(Node(this,nodeNum));

	numNodes++;
	assert(nodeNum >= 0);
	assert(nodeOrder[nodes[nodeNum].indexInOrder] == Node(this, nodeNum));
	return Node(this, nodeNum);
}


/**
 * Creates a new node with value one. If before there was no reduced node with value one, the new node is reduced.
 */
Node PowerCircuitGraph::newOneNode()
{
	Node node = newNode();
	if(numReducedNodes == 0)
	{
		moveNodeIntoReducedPart(node, 0);
	}
	return node;
}


/**
 * Creates a new marking with value one. If necessary a new node with value one is created.
 */
Marking PowerCircuitGraph::newOneMarking()
{
	if(numReducedNodes == 0)
		newOneNode();
	return newMarking(std::list< std::pair<Node,Sign> >(1, std::pair<Node,Sign>(nodeOrder[0],PLUS)));
}


/**
 * Creates a new marking with value zero.
 */
Marking PowerCircuitGraph::newZeroMarking()
{
	return newMarking();
}


/**
 * Creates a new marking with the same value as mark.
 */
Marking PowerCircuitGraph::newCopyMarking(const Marking& mark)
{
	assert(checkMarkingValid(mark));
	std::list< std::pair<Node,Sign> >& nodeList =markings[mark.id].nodes;
	return newMarking(nodeList);
}


/**
 * Creates a new marking consisting of zeros and only one one at position onePos (onePos indicates the position of the node
 * in the reduced order). OnePos should be in the range [0, numReducesNodes -1], otherwise the behavior is undefined.
 */
Marking PowerCircuitGraph::newUnitMarking(int onePos)
{
	return newMarking(std::list< std::pair<Node,Sign> >(1, std::pair<Node,Sign>(nodeOrder[onePos],PLUS)));
}


/**
 * Increases a list of signed nodes by one. The list has to be sorted (by index in nodeOrder) and all nodes
 * have to be in the reduced part of the circuit.
 */
std::list< std::pair<Node,Sign> > PowerCircuitGraph::inc(const std::list< std::pair<Node,Sign> >& nodeList)
{
	if(numReducedNodes==0)
		newOneNode();
	std::list< std::pair<Node,Sign> > result = nodeList;
	sortNodeList(result);

	std::list< std::pair<Node,Sign> >::iterator nodeIter = result.begin();
	unsigned int i = 0;
	while(nodeIter !=result.end() && nodeIter->second == PLUS && nodes[nodeIter->first.id].indexInOrder == i && nodes[nodeIter->first.id].BV)
	{
		 nodeIter = result.erase(nodeIter);
		 i++;
	}
	if(nodeIter == result.end() || nodes[nodeIter->first.id].indexInOrder != i)
	{
		assert(nodeOrder[i].id>=0);
		result.push_back(std::pair<Node,Sign>(nodeOrder[i],PLUS));
	}
	else if(nodeIter->second != PLUS)
	{
		if(nodeIter->second == ZERO)
			nodeIter->second = PLUS;
		else
			result.erase(nodeIter);
	}
	else
	{
		insertNewPowerOfTwoNode(i + 1);
		assert(nodeOrder[i+1].id>=0);
		result.push_back(std::pair<Node,Sign>(nodeOrder[i + 1], PLUS));
	}
	return result;
}


//-------------------------------------------------------------------------------------------
// Other auxiliary procedures dealing with reduced nodes/markings

/**
 * Moves a node to the reduced part of the matrix and marks it as reduced.
 * NOTICE: Its Lambda-marking has to be reduced already! Does NOT update the BV-vector!
 */
void PowerCircuitGraph::moveNodeIntoReducedPart(Node node, int newPos)
{
	unsigned int oldPos = nodes[node.id].indexInOrder;
	assert(oldPos >= numReducedNodes && newPos <= (int) numReducedNodes);

	for(int i = oldPos ; i > newPos ; i--)
	{
		nodeOrder[i] = nodeOrder[i - 1];
		nodes[nodeOrder[i].id].indexInOrder = i;
	}

	nodeOrder[newPos] = node;
	nodes[node.id].indexInOrder = newPos;
	assert(nodeOrder[nodes[node.id].indexInOrder] == node);
	numReducedNodes++;
}


/**
 * Compares two list of signed nodes. All nodes in the lists have to be in the reduced part of the circuit
 * if l1 < l2, -1 is returned
 * if l1 == l2, 0 is returned
 * if l1 > l2, +1 is returned
 * May have side effects.
 */
int PowerCircuitGraph::compare(std::list< std::pair<Node,Sign> >&  l1, std::list< std::pair<Node,Sign> >&  l2)
{
	sortNodeList(l1);
	sortNodeList(l2);

	std::list< std::pair<Node,Sign> > intersectionList;
	std::list< std::pair<Node,Sign> >::reverse_iterator iter1 = l1.rbegin();
	std::list< std::pair<Node,Sign> >::reverse_iterator iter2 = l2.rbegin();

	int val1;
	int carry = 0;
	int val2;
	int val=0;
	int oldNode=-1;
	int newNode;
	while(iter1 != l1.rend() || iter2 != l2.rend())
	{
		if(iter1 == l1.rend() || (iter2 != l2.rend() && compareNodesLessThan(*iter1,*iter2)) )
		{
			newNode=nodes[iter2->first.id].indexInOrder;
			val1=0;
			val2=signToInt(iter2->second);
			iter2++;
		}
		else if(iter2 == l2.rend() || compareNodesLessThan(*iter2,*iter1))
		{
			newNode=nodes[iter1->first.id].indexInOrder;
			val2=0;
			val1=signToInt(iter1->second);
			iter1++;
		}
		else
		{
			newNode=nodes[iter2->first.id].indexInOrder;
			val1=signToInt(iter1->second);
			val2=signToInt(iter2->second);

			iter1++;
			iter2++;
		}

		if(newNode==oldNode-1&&nodes[nodeOrder[newNode].id].BV)
			val= val1-val2+2*carry;
		else
			val= val1-val2+4*carry;

		if(val > 1)
			return 1;
		else if(val < -1)
			return -1;

		carry=val;
		oldNode=newNode;
	}
	return val;
}


/**
 * Sorts a list of signed nodes by there index in nodeOrder
 */
void PowerCircuitGraph::sortNodeList(std::list< std::pair<Node,Sign> >&  l)
{
	if(l.size()*4*log2(numNodes)<numNodes)
		l.sort(compareNodesLessThan);
	else
		l.sort(compareNodesLessThan); //ToDo to save a log-factor (to be theoretically more efficient)
}


//-------------------------------------------------------------------------------------------
// Auxiliary procedures for incMarking


/**
 * Creates a new node with value 2^power (hence with successor-marking of value power) and inserts it into the
 * reduced part of the PC. All nodes with values 2^i (0 <= i < power) have to exist already!!!
 */
void PowerCircuitGraph::insertNewPowerOfTwoNode(unsigned int power)
{
	if(nodes[nodeOrder[power-1].id].BV)
		return;
	Node node=newNode(calculateCompactRepresentation(power));

	moveNodeIntoReducedPart(node,power);
	nodes[nodeOrder[power-1].id].BV=true;

	if(power<numReducedNodes-1)
	{
		std::list< std::pair<Node,Sign> > compareList=calculateCompactRepresentation(power+1);
		if(compare(nodes[nodeOrder[power + 1].id].successors, compareList) == 0)
		{
			nodes[node.id].BV = true;
		}
		else
		{
			nodes[node.id].BV = false;
		}
	}
	else
	{
		nodes[node.id].BV = false;
	}
}


/**
 * Calculates the compact Representation of a given integer and returns it as signed list of nodes.
 */
std::list< std::pair<Node,Sign> > PowerCircuitGraph::calculateCompactRepresentation(int n)
{
	unsigned int temp=(n>=0)?n:(-n);
	int uebertrag=0;
	std::list< std::pair<Node,Sign> > result;
	for(int i=0;i<32;i++)
	{
		if(uebertrag+(temp&1)==0)
		{
			uebertrag=0;
		}
		else if(uebertrag+(temp&1)==2)
		{
			uebertrag=1;
		}
		else if((temp&2)!=0)
		{
			if(n>0)
				result.push_back(std::pair<Node,Sign>(nodeOrder[i],MINUS));
			else
				result.push_back(std::pair<Node,Sign>(nodeOrder[i],PLUS));
			uebertrag=1;
		}
		else
		{
			if(n>0)
				result.push_back(std::pair<Node,Sign>(nodeOrder[i],PLUS));
			else
				result.push_back(std::pair<Node,Sign>(nodeOrder[i],MINUS));
			uebertrag=0;
		}
		temp=temp>>1;
	}
	return result;
}


//------------------------------------------------------------------------------------------------
// Auxiliary procedures for reduce()


/**
 * Marks recursively all successor of nodes in nodeList and their successor in a vector marked.
 */
void PowerCircuitGraph::markSucessors(std::vector<bool> marked, std::list< std::pair<Node,Sign> >& nodeList)
{
	for(std::list< std::pair<Node,Sign> >::iterator nodeIter = nodeList.begin() ; nodeIter != nodeList.end() ; nodeIter++)
	{
		if(!marked[nodes[nodeIter->first.id].indexInOrder])
		{
			marked[nodes[nodeIter->first.id].indexInOrder] = true;
			markSucessors(marked, nodes[nodeIter->first.id].successors);
		}
	}
}


/**
 * Sets the BV-value of a node correctly. The node has to be in the reduced part of the PC.
 */
void PowerCircuitGraph::setBV(Node node)
{
	unsigned int nodePos=nodes[node.id].indexInOrder;
	if(nodePos<numReducedNodes-1)
	{
		std::list< std::pair<Node,Sign> > compareList = inc(nodes[node.id].successors);

		if(compare(compareList,nodes[nodeOrder[nodes[node.id].indexInOrder + 1].id].successors) == 0)
			nodes[node.id].BV = true;
		else
			nodes[node.id].BV = false;

	}
	else
	{
		nodes[node.id].BV=false;
	}
}


/**
 * Creates a new node with the double value of a reduced node "node" and inserts it into the reduced part of the power-circuit.
 */
Node PowerCircuitGraph::newDoubleNode(Node node)
{
	assert(nodeOrder[nodes[node.id].indexInOrder] == node);
	if(nodes[node.id].indexInOrder>=numReducedNodes)
	{
		print();
		std::cout<<"Error in newDoubleNode: node not reduced! RowIndex:"<<nodes[node.id].indexInOrder<<" numReducedNodes: "<<numReducedNodes<<std::endl;
		assert(false);
		return Node(this, -1);
	}
	if(nodes[node.id].BV) // if there is alreay a node with the double value
	{
		return nodeOrder[nodes[node.id].indexInOrder + 1];
	}
	Node newNode=this->newNode(inc(nodes[node.id].successors));

	moveNodeIntoReducedPart(newNode, nodes[node.id].indexInOrder + 1);
	nodes[node.id].BV=true;
	setBV(newNode);
	return newNode;
}


/**
 * Returns an iterator pointing on the element of nodeList, which refers to the node with index nodeIndex.
 * If there is no such node, nodeList.end() is returned.
 */
std::list< std::pair<Node,Sign> >::iterator PowerCircuitGraph::findNodeInList(unsigned int nodeIndex, std::list< std::pair<Node,Sign> >& nodeList)
{
	std::list< std::pair<Node,Sign> >::iterator nodeIter = nodeList.begin();
	while(nodeIter!= nodeList.end() && nodes[nodeIter->first.id].indexInOrder != nodeIndex)
		nodeIter++;
	return nodeIter;
}


/**
 * Checks all non-reduced columns if they use oldNode, and, if so, oldNode is replaced by the node newNode.
 * newNode has to be in the reduced part of the PC!!! If a marking uses both oldNode and newNode, new nodes are created
 * and inserted into the reduced part of the PC.
 */
void PowerCircuitGraph::removeDoubleNodesFromMarkings(Node oldNode, Node newNode, std::list<NodeUsedByType>* nodeUsedBy)
{
	Sign result;
	Sign carry;

	assert(oldNode.id != newNode.id);
	std::list<NodeUsedByType>::iterator newNodeIter = nodeUsedBy[newNode.id].begin();
	while(!nodeUsedBy[oldNode.id].empty())
	{
		int id = nodeUsedBy[oldNode.id].front().id;
		while(newNodeIter != nodeUsedBy[newNode.id].end() && newNodeIter->id < id )
			newNodeIter++;

		if(newNodeIter == nodeUsedBy[newNode.id].end() || newNodeIter->id > id)
		{
			nodeUsedBy[oldNode.id].front().position->first = newNode;
			nodeUsedBy[newNode.id].insert(newNodeIter, nodeUsedBy[oldNode.id].front());
		}
		else //newNodeIter->second == id
		{
			assert( newNodeIter->id == id );
			assert(newNodeIter->nodeList == nodeUsedBy[oldNode.id].front().nodeList);

			addSigns(newNodeIter->position->second, nodeUsedBy[oldNode.id].front().position->second, result, carry);

			assert(result == 0);

			if(carry == 0)
			{
				nodeUsedBy[oldNode.id].front().nodeList->erase(newNodeIter->position);
				newNodeIter = nodeUsedBy[newNode.id].erase(newNodeIter);

				nodeUsedBy[oldNode.id].front().nodeList->erase(nodeUsedBy[oldNode.id].front().position);
			}
			else //carry !=0
			{
				unsigned int i = nodes[newNodeIter->position->first.id].indexInOrder;
				nodeUsedBy[oldNode.id].front().nodeList->erase(newNodeIter->position);
				newNodeIter = nodeUsedBy[newNode.id].erase(newNodeIter);
				i++;
				while(true)
				{
					if(!nodes[nodeOrder[i - 1].id].BV)
					{
						Node node = newDoubleNode(nodeOrder[i - 1]);
						assert(node.id < (signed int)(3* nodes.size()));

						nodeUsedBy[oldNode.id].front().nodeList->push_front(std::pair<Node,Sign>(node, carry));
						nodeUsedBy[node.id].push_back(NodeUsedByType(nodeUsedBy[oldNode.id].front().nodeList, nodeUsedBy[oldNode.id].front().nodeList->begin(), id));

						break;
					}
					assert(i<numReducedNodes);
					std::list< std::pair<Node,Sign> >::iterator nodeIter = findNodeInList(i, *(nodeUsedBy[oldNode.id].front().nodeList));
					if(nodeIter == nodeUsedBy[oldNode.id].front().nodeList->end())
					{
						nodeUsedBy[oldNode.id].front().nodeList->push_front(std::pair<Node,Sign>(nodeOrder[i], carry));
						std::list<NodeUsedByType>::iterator tmpNodeIter = nodeUsedBy[nodeOrder[i].id].begin();
						while( tmpNodeIter != nodeUsedBy[nodeOrder[i].id].end() && tmpNodeIter->id < id )
							tmpNodeIter++;
						nodeUsedBy[nodeOrder[i].id].insert(tmpNodeIter, NodeUsedByType(nodeUsedBy[oldNode.id].front().nodeList, nodeUsedBy[oldNode.id].front().nodeList->begin(),id));
						break;
					}

					bool ready = (nodeIter->second != carry);
					assert(nodeIter->second!=ZERO);
					nodeUsedBy[nodeOrder[i].id].remove(NodeUsedByType(nodeUsedBy[oldNode.id].front().nodeList, nodeIter, id));
					nodeUsedBy[oldNode.id].front().nodeList->erase(nodeIter);

					if(ready)
					{
						break;
					}
					i++;
				}
				nodeUsedBy[oldNode.id].front().nodeList->erase(nodeUsedBy[oldNode.id].front().position);
			}
		}
		nodeUsedBy[oldNode.id].pop_front();
	}
}


/**
 * Returns the position, where a given node should be inserted in the reduced part of a circuit.
 * If there is already a reduced node with the same value, equal is set to true, otherwise it is set to false.
 */
int PowerCircuitGraph::findNewPosOfReducedNode(Node node, bool& equal)
{
	int first=0;
	int last=numReducedNodes;
	int next;
	if(numReducedNodes==0)
		return 0;

	//now we have last - first >= 1 this remains true through the whole procedure
	std::list< std::pair<Node,Sign> >& compareList = nodes[node.id].successors;
	int compResult=compare(compareList, nodes[nodeOrder[first].id].successors);

	equal=false;
	if(compResult<0)
		return 0;
	else if(compResult==0)
	{
		equal = true;
		return 0;
	}
	next = (first + last) / 2;

	//now the value of col ist always bigger than the value of first and smaller than the value of last
	while(first!=next)//that means last - first > 1
	{
		compResult=compare(compareList,nodes[nodeOrder[next].id].successors);
		if(compResult==0)
		{
			equal=true;
			return next;
		}
		else if(compResult<0)
		{
			last = next;
		}
		else
		{
			first = next;
		}
		next = (first + last) / 2;
	}
	return last;
}


/**
 * Inserts the node into the reduced part of the PC.
 * Its successor-marking has to have support in reduced part of the PC!
 */
int PowerCircuitGraph::insertNodeIntoReduced(Node node, std::list<NodeUsedByType>* nodeUsedBy)
{
	bool equal;
	int newPos = findNewPosOfReducedNode(node, equal);
	assert(nodeOrder[nodes[node.id].indexInOrder] == node);
	assert(nodes[node.id].indexInOrder >= numReducedNodes);
	assert(newPos <= (int)numReducedNodes);

	if(!equal)// if there is no reduced node with the same Lambda-Marking as node
	{
		moveNodeIntoReducedPart(node, newPos);
		if(newPos>0)
			setBV(nodeOrder[newPos-1]);
		setBV(node);
	}
	else // if there is already a reduced node with the same Lambda-Marking
	{
		removeDoubleNodesFromMarkings(node, nodeOrder[newPos], nodeUsedBy);
		deleteNode(node);
	}
	return newPos;
}


/**
 * Recursive procedure to sort topologically starting with node node.
 * This procedure has to be called for all nodes which should be sorted topologically.
 */
void PowerCircuitGraph::topSortNode(Node node, std::vector<bool>& visited, std::vector<Node>& topSortPerm)
{
	visited[nodes[node.id].indexInOrder]=true;
	std::list< std::pair<Node,Sign> >& nodeList = nodes[node.id].successors;

	for(std::list< std::pair<Node,Sign> >::iterator iter = nodeList.begin(); iter!= nodeList.end() ; iter++)
	{
		if(!visited[nodes[iter->first.id].indexInOrder])
		{
			assert(iter->second!=ZERO);
			topSortNode(iter->first, visited, topSortPerm);
		}

	}
	topSortPerm.push_back(node);
}



/**
 *  Moves all nodes in nodeList and markings in markingList to the reduced part of the PC.
 *  At the end nodeList is the empty list!!!
 *  Nodes in nodeList MAY NOT be reduced already!!
 *  The support of markings in markingList and successor-markings of node in nodeList has to be in nodeList
 *  or the already reduced part of the PC.
 */
void PowerCircuitGraph::extendTree(std::vector<Node>& nodeList, std::vector<Marking>& markingList)
{
	std::vector<Node>  topSortPerm;
	int numNodesToSort = nodeList.size();
	topSortPerm.reserve(numNodesToSort);

	std::list<NodeUsedByType> * nodeUsedBy = new std::list<NodeUsedByType> [3 * nodes.size()];

	//reserve enough storage capacity, to assure that during the whole extendTree-procedure no reallocation is required.
	nodes.reserve(std::max(3 * numNodes +1, (unsigned int)nodes.capacity()));

	std::vector<bool>  visited;
	visited.assign(numNodes, true);
	int id = 0;

	//initialize nodeUsedBy (pointers from nodes to markings)
	for(std::vector<Node>::iterator nodeIter = nodeList.begin() ; nodeIter != nodeList.end() ; nodeIter++)
	{
		visited[nodes[nodeIter->id].indexInOrder] = false;

		std::list< std::pair<Node,Sign> >& nodeList2 = nodes[nodeIter->id].successors;
		for(std::list< std::pair<Node,Sign> >::iterator nodeIter2 = nodeList2.begin() ; nodeIter2 != nodeList2.end() ; nodeIter2++)
		{
			nodeUsedBy[nodeIter2->first.id].push_back(NodeUsedByType(&nodeList2,nodeIter2,id));
		}
		id++;
	}

	for(std::vector<Marking>::iterator markIter = markingList.begin() ;  markIter != markingList.end() ; markIter++)
	{
		std::list< std::pair<Node,Sign> >& nodeList2 = markings[markIter->id].nodes;

		for(std::list< std::pair<Node,Sign> >::iterator nodeIter2 = nodeList2.begin() ; nodeIter2 != nodeList2.end() ; nodeIter2++)
		{
			nodeUsedBy[nodeIter2->first.id].push_back(NodeUsedByType(&nodeList2,nodeIter2,id));
		}
		id++;
	}

	//sort topologically
	for(std::vector<Node>::iterator nodeIter = nodeList.begin() ; nodeIter != nodeList.end() ; nodeIter++)
	{
		if((!visited[nodes[nodeIter->id].indexInOrder]))
			topSortNode(*nodeIter, visited, topSortPerm);
	}
	assert((int)topSortPerm.size()==numNodesToSort);

	for(int i=0; i < numNodesToSort ; i++)
	{
		insertNodeIntoReduced(topSortPerm[i], nodeUsedBy);
	}

	for(std::vector<Marking>::iterator markIter = markingList.begin() ;  markIter != markingList.end() ; markIter++)
		markings[markIter->id].reduced = true;
	delete[] nodeUsedBy;
}




//-------------------------------------------------------------------------------------------------
// implementation of private abstract methods of PowerCircuit


/**
 * Increases the reference counter of a a marking. (Called only by constructors and assignment-operator of markings)
 */
void PowerCircuitGraph::incMarkingRefCount(const Marking& mark)
{
	assert(checkMarkingValid(mark));
	markings[mark.id].refCount++;
}

/**
 * Decreases the reference counter of a a marking, and if required deletes a marking.
 * (Called only by destructor and assignment-operator of markings)
 */
void PowerCircuitGraph::decMarkingRefCount(Marking& mark)
{
	assert(checkMarkingValid(mark));
	if(markings[mark.id].refCount > 1)
		markings[mark.id].refCount--;
	else
		deleteMarking(mark);
}


/**
 * Creates a new Marking with the inverse value of mark.
 */
Marking PowerCircuitGraph::invMarking(const Marking& mark)
{
	assert(checkMarkingValid(mark));
	std::list< std::pair<Node,Sign> > invList;
	for(std::list< std::pair<Node,Sign> >::iterator nodeIter = markings[mark.id].nodes.begin()  ; nodeIter != markings[mark.id].nodes.end() ; nodeIter++)
	{
		if(nodeIter->second == PLUS)
			invList.push_back(std::pair<Node,Sign>(nodeIter->first,MINUS));
		else if(nodeIter->second == MINUS)
			invList.push_back(std::pair<Node,Sign>(nodeIter->first,PLUS));
	}
	return newMarking(invList);
}


/**
 * Creates a new Marking with the sum of mark1 and mark2. The user has to assure that for each node,
 * the sum lies in the range [-1,+1]!!!
 */
Marking PowerCircuitGraph::addMarkings(const Marking& mark1, const Marking&  mark2)
{

	checkMarkingValid(mark1);
	checkMarkingValid(mark2);
	std::list< std::pair<Node,Sign> >&  nodeList1 = markings[mark1.id].nodes;
	std::list< std::pair<Node,Sign> >&  nodeList2 = markings[mark2.id].nodes;
	sortNodeList(nodeList1);
	sortNodeList(nodeList2);

	std::list< std::pair<Node,Sign> > addList;
	std::list< std::pair<Node,Sign> >::iterator nodeIter1 = nodeList1.begin();
	std::list< std::pair<Node,Sign> >::iterator nodeIter2 = nodeList2.begin();
	while(nodeIter1 != nodeList1.end() || nodeIter2 != nodeList2.end())
	{
		if(nodeIter2 == nodeList2.end() || (nodeIter1 != nodeList1.end() && nodes[nodeIter1->first.id].indexInOrder < nodes[nodeIter2->first.id].indexInOrder))
		{
			addList.push_back(*nodeIter1);
			nodeIter1++;
		}
		else if(nodeIter1 == nodeList1.end() || nodes[nodeIter1->first.id].indexInOrder > nodes[nodeIter2->first.id].indexInOrder)
		{
			addList.push_back(*nodeIter2);
			nodeIter2++;
		}
		else
		{
			Sign sum, carry;
			addSigns(nodeIter1->second,nodeIter2->second,sum,carry);
			assert(carry==ZERO);
			addList.push_back(std::pair<Node,Sign>(nodeIter1->first,sum));
			nodeIter1++;
			nodeIter2++;
		}
	}
	return newMarking(addList);
}


/**
 * Creates a new Marking m with the intersection of mark1 and mark2, that is
 * m(p) = mark1(p), 	if (mark1(p)==mark2(p))
 * m(p) = 0 			else
 * for every node p.
 */
Marking PowerCircuitGraph::intersectMarkings(const Marking& mark1, const Marking&  mark2)
{
	checkMarkingValid(mark1);
	checkMarkingValid(mark2);
	std::list< std::pair<Node,Sign> >&  nodeList1 = markings[mark1.id].nodes;
	std::list< std::pair<Node,Sign> >&  nodeList2 = markings[mark2.id].nodes;
	sortNodeList(nodeList1);
	sortNodeList(nodeList2);

	std::list< std::pair<Node,Sign> > intersectionList;
	std::list< std::pair<Node,Sign> >::iterator nodeIter1 = nodeList1.begin();
	std::list< std::pair<Node,Sign> >::iterator nodeIter2 = nodeList2.begin();
	while(nodeIter1 != nodeList1.end() && nodeIter2 != nodeList2.end())
	{
		if(nodeIter1->first == nodeIter2->first)
			intersectionList.push_back(*nodeIter1);
		nodeIter1++;
		nodeIter2++;
	}

	return newMarking(intersectionList);
}


/**
 * Increases a given marking by one. The result is NOT reduced. If mark is reduced, the result has support in the reduced
 * part of the circuit.
 */
Marking PowerCircuitGraph::incMarking(const Marking& mark)
{
	assert(checkMarkingValid(mark));

	if(!markings[mark.id].reduced)
	{
		Marking result=newCopyMarking(mark);
		Node oneNode = newOneNode();
		markings[result.id].nodes.push_back(std::pair<Node,Sign>(oneNode,PLUS));
		return result;
	}
	else
	{
		return newMarking(inc(markings[mark.id].nodes));
	}
}


/**
 * Clones a marking. That means all nodes in the support of the marking are cloned and a new marking with
 * the same value as the original one (using only the newly created nodes) is created.
 */
Marking PowerCircuitGraph::cloneMarking(const Marking& mark)
{
	assert(checkMarkingValid(mark));
	Marking result=newZeroMarking();

	for(std::list< std::pair<Node,Sign> >::iterator nodeIter = markings[mark.id].nodes.begin()  ; nodeIter != markings[mark.id].nodes.end() ; nodeIter++)
	{
		assert(nodeIter->second!=ZERO);
		Node node=cloneNode(nodeIter->first);
		markings[result.id].nodes.push_back( std::pair<Node,Sign> (node, nodeIter->second));

	}
	return result;
}


/**
 * Clones a node, that means a new node is inserted with a copy of the successor-marking.
 * Even if the old node was reduced, the new one is not!
 */
Node PowerCircuitGraph::cloneNode(Node node)
{
	return newNode(nodes[node.id].successors);
}


/**
 * Return a deep copy of a marking. (that means a new marking with the same nodes as the old one)
 */
Marking	PowerCircuitGraph::copyMarking(const Marking& m)
{
	assert(checkMarkingValid(m));
	return newCopyMarking(m);
}


/**
 * Returns whether a marking is reduced.
 */
bool PowerCircuitGraph::isMarkingReduced(const Marking& m) const
{
	checkMarkingValid( m);
	return markings[m.id].reduced;
}


/**
 * Compares the values of two markings. Returns
 * -1, if m1<m2,
 * 0,  if m1==m2
 * 1,  if m1>m1
 * Only reduced markings may be compared.
 */
int	PowerCircuitGraph::compareMarkings(const Marking& m1, const Marking&  m2)
{
	assert(checkMarkingValid(m1));
	assert(checkMarkingValid(m2));
	if(!markings[m1.id].reduced||!markings[m2.id].reduced)
	{
		throw "Only markings in the reduced part of the circuit can be compared.";
	}
	return compare(markings[m1.id].nodes, markings[m2.id].nodes);
}


/**
 *	Returns the index of a node in the reduces part of the PC
 *	(thus the number of smaller nodes in the reduced part).
 *	If the node is not in the reduced part, -1 is returned.
 */
int	PowerCircuitGraph::getRedNodeOrd(Node n)
{
	assert(checkNodeValid(n));
	if(nodes[n.id].indexInOrder<numReducedNodes)
		return nodes[n.id].indexInOrder;
	else
		return -1;
}


/**
 * Returns the sign of a given node n in a marking m. 1 stands for PLUS, -1 for MINUS.
 */
int	PowerCircuitGraph::getNodeSignInMarking(Node n, const Marking&  m) const
{
	assert(checkNodeValid(n));
	assert(checkMarkingValid(m));

	for(std::list< std::pair<Node,Sign> >::const_iterator nodeIter = markings[m.id].nodes.begin()  ; nodeIter != markings[m.id].nodes.end() ; nodeIter++)
	{
		assert(nodeIter->second!=ZERO);
		if(nodeIter->first.id == n.id)
			return compareSigns(nodeIter->second,ZERO);
	}
	return 0;
}


/**
 * Returns whether the marking m is a successor-marking.
 */
bool PowerCircuitGraph::isSuccessorMarking(const Marking& m) const
{
	assert(checkMarkingValid(m));
	return false;
}


/**
 * Takes as input a marking m. If m is a successor-marking,
 * the node whose successor-marking it is is returned.
 * Otherwise the undefined node is returned.
 */
Node PowerCircuitGraph::getIncidentNode(const Marking& m)
{
	assert(checkMarkingValid(m));
	return Node(this,-1);
}


/**
 * Returns the smallest non-zero node of a reduced marking m.
 * If m is the zero-marking the undefined node is returned.
 * If m is not reduced the it might return another node than the smallest.
 */
Node PowerCircuitGraph::getSmallestNodeInMarking(const Marking& m)
{
	assert(checkMarkingValid(m));
	Node minNode;
	unsigned int min = numNodes;

	for(std::list< std::pair<Node,Sign> >::iterator nodeIter = markings[m.id].nodes.begin()  ; nodeIter != markings[m.id].nodes.end() ; nodeIter++)
	{
		assert(nodeIter->second!=ZERO);
		assert(nodeOrder[nodes[nodeIter->first.id].indexInOrder] == nodeIter->first);
		assert(nodes[nodeIter->first.id].indexInOrder<numNodes);
		assert(nodes[nodeIter->first.id].indexInOrder >= 0);
		if(nodes[nodeIter->first.id].indexInOrder < min)
		{
			min = nodes[nodeIter->first.id].indexInOrder;
			minNode = nodeIter->first;
		}
	}
	if(min != numNodes)
		return minNode;
	else
		return Node(this,-1);
}


/**
 * Returns a copy of the successor-marking of the node n. (thus the returned marking is not a successormarking)
 */
Marking	PowerCircuitGraph::getSuccMarking(Node n)
{
	assert(checkNodeValid(n));
	return newMarking(nodes[n.id].successors);
}


//---------------------------------------------------------------------------------------------
// public methods


/**
 * Reduces the whole PC and all its markings.
 */
void PowerCircuitGraph::reduce()
{
	std::vector<Node>  nodeList;
	nodeList.reserve(numNodes-numReducedNodes);
	std::vector<Marking> markingList;
	markingList.reserve(numMarkings);

	//put nodes, which should be reduced (i.e. all non-reduced nodes) into temporary list
	for(unsigned int i = numReducedNodes ; i < numNodes ; i++)
	{
		nodeList.push_back(nodeOrder[i]);
	}

	for(unsigned int j = 0 ; j<markings.size() ; j++)
	{
		if(!markings[j].deleted && !markings[j].reduced)
			markingList.push_back(Marking(this,j));
	}

	extendTree(nodeList, markingList);
}


/**
 * Reduces all the markings in markingList, and the therefore required nodes.
 */
void PowerCircuitGraph::reduce(std::list<Marking> markingList)
{
	std::vector<Node>  nodeVector;
	nodeVector.reserve(numNodes-numReducedNodes);
	std::vector<Marking> markingVector;
	std::vector<bool> needed;
	needed.assign(numNodes-numReducedNodes, false);
	markingVector.reserve(markingList.size());

	for(std::list<Marking>::iterator markIter = markingList.begin() ; markIter != markingList.end(); markIter++)
	{
		assert(checkMarkingValid(*markIter));
		markSucessors(needed, markings[markIter->id].nodes);
		if(!markings[markIter->id].reduced)
			markingVector.push_back(*markIter);
	}

	//put nodes, which should be reduced (i.e. all non-reduced nodes) into temporary list
	for(unsigned int i = numReducedNodes ; i < numNodes ; i++)
	{
		if(needed[i])
			nodeVector.push_back(nodeOrder[i]);
	}

	extendTree(nodeVector, markingVector);
}


/**
 *  Moves all nodes in nodeList and markings in markingList to the reduced part of the PC.
 *  Nodes in nodeList MAY NOT be reduced already!!
 *  The support of markings in markingList has to be in nodeList or the already reduced part of the PC.
 */
void PowerCircuitGraph::reduce(std::vector<Node>& nodeVector, std::vector<Marking>& markingVector)
{
	for(std::vector<Node>::iterator nodeIter=nodeVector.begin() ; nodeIter!=nodeVector.end(); nodeIter++)
	{
		assert(nodes[nodeIter->id].indexInOrder >= numReducedNodes);
	}
	for(std::vector<Marking>::iterator markIter=markingVector.begin() ;markIter!=markingVector.end(); markIter++)
	{
		assert(!markings[markIter->id].reduced);
	}

	extendTree(nodeVector, markingVector);
}


/**
 * Creates a new marking with value val. The therefore required nodes are inserted into the reduced part of the PC.
 */
Marking	PowerCircuitGraph::createMarking(int val)
{
	if(val==0)
		return newZeroMarking();
	else if(val==1)
		return newOneMarking();
	else
	{
		Marking m =newZeroMarking();
		int temp;
		Sign sign;
		if(val<0)
		{
			sign=MINUS;
			temp=-val;
		}
		else
		{
			sign=PLUS;
			temp=val;
		}

		int i=0;
		if(numReducedNodes==0)
			newOneNode();
		if((temp&0x1)==1)
			markings[m.id].nodes.push_back(std::pair<Node,Sign>(nodeOrder[0], sign));
		temp=temp>>1;
		i++;

		while(temp!=0)
		{
			assert((unsigned int)i<8*sizeof(int));
			if(!nodes[nodeOrder[i-1].id].BV)
				insertNewPowerOfTwoNode(i);

			if((temp&0x1)==1)
				markings[m.id].nodes.push_back(std::pair<Node,Sign>(nodeOrder[i], sign));
			temp=temp>>1;
			i++;
		}
		return m;
	}
}


/**
 *  Creates a new marking, with the nodes in nodeList set to 1.
 */
Marking	PowerCircuitGraph::createMarking(const std::list<Node>& nodeList)
{
	Marking mark=newZeroMarking();
	for(std::list<Node>::const_iterator nodeIter = nodeList.begin() ; nodeIter != nodeList.end() ; nodeIter++)
	{
		assert(checkNodeValid(*nodeIter));
		markings[mark.id].nodes.push_back(std::pair<Node,Sign>(*nodeIter,PLUS));
	}
	return mark;
}


/**
 *  Creates a new marking, with the nodes in nodeList set to the respective sign.
 */
Marking	PowerCircuitGraph::createMarking(const std::list< std::pair<Node,Sign> >& nodeList )
{
	return newMarking(nodeList);
}


/**
 *  Creates a new marking, with the nodes in nodeList set to the respective sign in signList.
 */
Marking	PowerCircuitGraph::createMarking(const std::list<Node>& nodeList, const std::list<int>& signList)
{
	Marking mark=newZeroMarking();
	std::list<int>::const_iterator signIter = signList.begin();
	for(std::list<Node>::const_iterator nodeIter = nodeList.begin() ;	nodeIter != nodeList.end() &&  signIter != signList.end(); nodeIter++)
	{
		if(*signIter > 0)
			markings[mark.id].nodes.push_back(std::pair<Node,Sign>(*nodeIter,PLUS));
		else
			markings[mark.id].nodes.push_back(std::pair<Node,Sign>(*nodeIter,MINUS));

		signIter++;
	}
	return mark;
}


/**
 * Creates a new node with a copy of m as successor-marking.
 */
Node PowerCircuitGraph::createNode(const Marking& m)
{
	assert(checkMarkingValid(m));
	return newNode(markings[m.id].nodes);
}


/**
 * Adds the marking p to all successor-markings of nodes in the support of m. (So, if supp(m) and supp(p) are
 * disjoint, that means, that arrows are drawn from all nodes in supp(m) to p).
 * The user has to assure that, for each node, the sum lies in the range [-1,+1] (for example using clone())!!!
 * The user is required to check that the resulting graph has no cycles and no successor-marking becomes negative.
 * All nodes in supp(m) have to be non-reduced (use clone() to assure that).
 */
void PowerCircuitGraph::connect(const Marking& m, const Marking&  p)
{
	assert(checkMarkingValid(m));
	checkMarkingValid(p);

	std::list< std::pair<Node,Sign> >&  nodeListM = markings[m.id].nodes;
	std::list< std::pair<Node,Sign> >&  nodeListP = markings[p.id].nodes;
	sortNodeList(nodeListP);

	for(std::list< std::pair<Node,Sign> >::iterator nodeIterM = nodeListM.begin() ; nodeIterM != nodeListM.end() ; nodeIterM++ )
	{
		if(nodes[nodeIterM->first.id].indexInOrder<numReducedNodes)
		{
			throw "Can only connect non-reduced nodes.";
		}
		std::list< std::pair<Node,Sign> >&  nodeListTarget = nodes[nodeIterM->first.id].successors;
		sortNodeList(nodeListTarget);

		std::list< std::pair<Node,Sign> >::iterator nodeIterTarget = nodeListTarget.begin();
		std::list< std::pair<Node,Sign> >::iterator nodeIterP = nodeListP.begin();
		while(nodeIterTarget != nodeListTarget.end() && nodeIterP != nodeListP.end())
		{
			if(nodes[nodeIterTarget->first.id].indexInOrder < nodes[nodeIterP->first.id].indexInOrder)
			{
				nodeIterTarget++;
			}
			else if(nodes[nodeIterTarget->first.id].indexInOrder > nodes[nodeIterP->first.id].indexInOrder)
			{
				nodeListP.insert(nodeIterTarget,*nodeIterP);
				nodeIterP++;
			}
			else
			{
				Sign sum, carry;
				addSigns(nodeIterTarget->second,nodeIterP->second,sum,carry);
				assert(carry==ZERO);
				if(sum==ZERO)
				{
					nodeIterTarget = nodeListTarget.erase(nodeIterTarget);
				}
				else
				{
					nodeIterTarget->second = sum;
					nodeIterTarget++;
				}
				nodeIterP++;
			}
		}
		if(nodeIterP != nodeListP.end())
		{
			assert(nodeIterTarget == nodeListTarget.end());
			nodeListTarget.insert(nodeIterTarget,nodeIterP,nodeListP.end());
		}


	}
}


/**
 * The same as connect, only that the marking p is subtracted instead of added.
 */
void PowerCircuitGraph::connectInv(const Marking& m, const Marking&  p)
{
	assert(checkMarkingValid(m));
	checkMarkingValid(p);


	std::list< std::pair<Node,Sign> >&  nodeListM = markings[m.id].nodes;
	std::list< std::pair<Node,Sign> >&  nodeListP = markings[p.id].nodes;
	sortNodeList(nodeListP);

	for(std::list< std::pair<Node,Sign> >::iterator nodeIterM = nodeListM.begin() ; nodeIterM != nodeListM.end() ; nodeIterM++ )
	{
		if(nodes[nodeIterM->first.id].indexInOrder<numReducedNodes)
		{
			throw "Can only connect non-reduced nodes.";
		}
		std::list< std::pair<Node,Sign> >&  nodeListTarget = nodes[nodeIterM->first.id].successors;
		sortNodeList(nodeListTarget);

		std::list< std::pair<Node,Sign> >::iterator nodeIterTarget = nodeListTarget.begin();
		std::list< std::pair<Node,Sign> >::iterator nodeIterP = nodeListP.begin();
		while(nodeIterTarget != nodeListTarget.end() && nodeIterP != nodeListP.end())
		{
			if(nodes[nodeIterTarget->first.id].indexInOrder < nodes[nodeIterP->first.id].indexInOrder)
			{
				nodeIterTarget++;
			}
			else if(nodes[nodeIterTarget->first.id].indexInOrder > nodes[nodeIterP->first.id].indexInOrder)
			{
				nodeListP.insert(nodeIterTarget,std::pair<Node,Sign>(nodeIterP->first,negateSign(nodeIterP->second)));
				nodeIterP++;
			}
			else
			{
				Sign sum, carry;
				addSigns(nodeIterTarget->second,negateSign(nodeIterP->second),sum,carry);
				assert(carry==ZERO);
				if(sum==ZERO)
				{
					nodeIterTarget = nodeListTarget.erase(nodeIterTarget);
				}
				else
				{
					nodeIterTarget->second = sum;
					nodeIterTarget++;
				}
				nodeIterP++;
			}
		}
		if(nodeIterP != nodeListP.end())
			assert(nodeIterTarget == nodeListTarget.end());
		while(nodeIterP != nodeListP.end())
		{
			nodeListTarget.push_back(std::pair<Node,Sign>(nodeIterP->first,negateSign(nodeIterP->second)));
			nodeIterP++;
		}
	}
}


/**
 * Returns the node with index ord (that means there are ord smaller nodes in the reduced part of the PC)
 * If there are less than ord nodes, the undefined node is returned.
 */
Node PowerCircuitGraph::getReducedNode(unsigned int ord)
{
	if(ord<numReducedNodes)
		return nodeOrder[ord];
	else
		return Node(this,-1);
}


/**
 * Returns a list of all nodes in the PowerCircuit.
 */
std::list<Node>	PowerCircuitGraph::getNodes()
{
	std::list<Node> nodes;
	for(unsigned int i=0;i<numNodes;i++)
	{
		nodes.push_back(nodeOrder[i]);
	}
	return nodes;
}


/**
 * Returns a list of all markings in the PowerCircuit.
 */
std::list<Marking>	PowerCircuitGraph::getMarkings()
{
	std::list<Marking> marks;
	for(unsigned int j=0;j<markings.size();j++)
	{
		if(!markings[j].deleted)
			marks.push_back(Marking(this, j));
	}
	return marks;
}


/**
 * Returns a list of all nodes which are non-zero in the marking m.
 */
std::list<Node>	PowerCircuitGraph::getMarkingNodes(const Marking& m)
{
	assert(checkMarkingValid(m));
	std::list<Node> nodes;
	for(std::list< std::pair<Node,Sign> >::iterator nodeIter = markings[m.id].nodes.begin()  ; nodeIter != markings[m.id].nodes.end() ; nodeIter++)
	{
		assert(nodeIter->second!=ZERO);
		nodes.push_back(nodeIter->first);
	}
	return nodes;
}


/**
 * Prints the adjacency-matrix of the underlying graph.
 */
void PowerCircuitGraph::print(std::ostream& os)
{
	os<<"numMarkings: "<<numMarkings<<std::endl;
	os<<"numNodes: "<<numNodes<<std::endl;

	std::vector<std::list<std::pair<Node,Sign> >::iterator > iterVector;
	for(unsigned int i = 0; i < numNodes; i++)
	{
		iterVector.push_back(nodes[nodeOrder[i].id].successors.begin());
	}
	for(unsigned int i = 0; i < numNodes; i++)
	{
		if(i<numReducedNodes)
			os << std::setw(3) << "t"<< (nodes[nodeOrder[i].id].BV?1:0);
		else
			os << std::setw(4) << " ";
		for (unsigned int j = 0 ; j < numNodes; j++)
		if(iterVector[j] != nodes[nodeOrder[j].id].successors.end() && iterVector[j]->first == nodeOrder[i])
		{
			os << std::setw(4) << signToInt(iterVector[j]->second);
			iterVector[j]++;
		}
		else
		{
			os << std::setw(4) << "0";
		}
		os << std::endl;
	}

	for(unsigned int i = 0 ; i < numMarkings ; i++)
	{
		if(!markings[i].deleted)
		{
			os << "M" << i << ": ";
			printMarking(markings[i].nodes);
			std::cout<<std::endl;
		}
	}
	os << std::endl;
	os.flush();
}


/**
 * Prints the nodes in a list.
 */
void PowerCircuitGraph::printMarking(const std::list< std::pair<Node, Sign> >& l ,std::ostream& os)
{
	for(std::list< std::pair<Node, Sign> >::const_iterator nodeIter = l.begin() ; nodeIter != l.end() ; nodeIter++)
	{
		os<<nodes[nodeIter->first.id].indexInOrder<<": "<<signToInt(nodeIter->second)<<", ";
	}
}


/**
 * Prints the nodes belonging to a marking.
 */
void PowerCircuitGraph::printMarking(const Marking& m,std::ostream& os)
{
	for(std::list< std::pair<Node, Sign> >::iterator nodeIter = markings[m.id].nodes.begin() ; nodeIter != markings[m.id].nodes.end() ; nodeIter++)
	{
		os<<"Node "<<nodes[nodeIter->first.id].indexInOrder<<" (id=="<<nodeIter->first.id<<"): "<<signToInt(nodeIter->second)<<std::endl;
	}
}


/**
 * Deletes all nodes of the marking. (So if other markings use those
 * nodes, their values may be changed!!!)
 */
 void PowerCircuitGraph::remove(const Marking& m)
 {
	 assert(checkMarkingValid(m));
	 std::list<Node> nodeList=getMarkingNodes(m);
	 while(!nodeList.empty())
	 {
		 Node node = nodeList.front();
		 nodeList.pop_front();
		 deleteNode(node);
	 }
 }


 /**
  * Makes a copy of a PowerCircuit. Only the markings in the markingsToClone list are copied.
  * Reduced nodes and markings are kept reduced.
  * Has to be deallocated via delete, once it is not needed anymore!!!
  */
 PowerCircuit* PowerCircuitGraph::clone(std::vector<Marking>& markingsToKeep)
 {
	PowerCircuitGraph* result = new PowerCircuitGraph();

	result->numReducedNodes=numReducedNodes;
	result->numNodes = numNodes;
	result->nodes = nodes;
	result->nodeOrder.resize(nodeOrder.size());
	for(unsigned int i = 0 ; i < numNodes; i++)
	{
		result->nodeOrder[i] = Node(result, nodeOrder[i].id);
		result->nodes[nodeOrder[i].id].successors.clear();
		for(std::list< std::pair<Node, Sign> >::iterator nodeIter = nodes[nodeOrder[i].id].successors.begin() ; nodeIter != nodes[nodeOrder[i].id].successors.end() ; nodeIter++)
			result->nodes[nodeOrder[i].id].successors.push_back(std::pair<Node, Sign>(Node(result, nodeIter->first.id), nodeIter->second));
	}
	result->firstDeletedNode = firstDeletedNode;

	//allocate the copied markings
	for(std::vector<Marking>::iterator markIter=markingsToKeep.begin() ; markIter != markingsToKeep.end(); markIter++)
	{
		std::list< std::pair<Node, Sign> > newNodeList;
		for(std::list< std::pair<Node, Sign> >::iterator nodeIter = markings[markIter->id].nodes.begin() ; nodeIter != markings[markIter->id].nodes.end() ; nodeIter++)
			newNodeList.push_back(std::pair<Node, Sign>(Node(result,nodeIter->first.id),nodeIter->second));
		Marking m=result->newMarking(newNodeList);
		*markIter=m;
	}
	return result;
 }


/**
 * Prints some statistical data.
 */
void PowerCircuitGraph::printStatistics(std::ostream& os)
{
	os<<"In this PC there are "<<getNumMarkings()<<" markings and "<<getNumNodes()<<" nodes."<<std::endl;
	os<<"Edges in underlying graph: "<<getNumEdges()<<". "<<std::endl<<std::endl;
}



bool PowerCircuitGraph::checkCyclesRecursive(int n, std::vector<bool>&  visited)
{
	if(visited[n])
	{
		draw("Cyclegraph");
		return true;
	}
	visited[n] = true;
	for(std::list< std::pair<Node,Sign> >::iterator nodeIter = nodes[nodeOrder[n].id].successors.begin() ;  nodeIter != nodes[nodeOrder[n].id].successors.end() ; nodeIter++)
	{
		assert(nodeIter->second != ZERO);
		if(checkCyclesRecursive(nodes[nodeIter->first.id].indexInOrder, visited))
			return true;

	}
	visited[n] = false;
	return false;
}

/**
 * Checks if the power-circuit has cycles.
 */
bool PowerCircuitGraph::checkCycles()
{
	std::vector<bool> visited;
	visited.assign(numNodes,false);
	for(unsigned int i = 0 ; i < numNodes ; i++)
	{
		if(checkCyclesRecursive(i, visited))
			return true;
	}
	return false;
}


/**
 * Checks if the power-circuit is consistent (i.e. if there are cycles).
 */
bool PowerCircuitGraph::checkConsistency()
{
	std::cout<<"checkConsistency"<<std::endl;
	if (checkCycles())
	{
		std::cout<<"checkConsistency: cycle"<<std::endl;
		print();
		return false;
	}
	for(int i = 0 ; i < (int)numReducedNodes - 1; i++)
	{
		if(compare(nodes[nodeOrder[i].id].successors, nodes[nodeOrder[i + 1].id].successors) != -1)
		{
			std::cout<<"checkConsistency: wrong order of nodes"<<std::endl;
			print();
			return false;
		}
	}

	if(getReducedNode(0).isDefined())
	{
		std::vector<Marking> markVector;
		std::list<Marking> markList = getMarkings();
		markVector.insert(markVector.end(), markList.begin() , markList.end());
		PowerCircuit* pcClone = this->clone(markVector);
		Marking m0 = pcClone->createMarking(0);
		Marking mNode = pcClone->getReducedNode(0).getSuccessorMarking();
		pcClone->reduce();

		if(mNode < m0)
		{
			std::cout<<"checkConsistency: negative successor-marking"<<std::endl;
			print();
			return false;
		}
	}
	return true;
}

}
