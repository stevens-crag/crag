/*
 * PowerCircuitCompMatrix.cpp
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
#include "SignMatrix.h"
#include "PowerCircuit.h"
#include "PowerCircuitCompMatrix.h"



namespace PC
{

// Initialize static members of PowerCircuitCompMatrix
PowerCircuitCompMatrix::MarkingVector 	PowerCircuitCompMatrix::markings;
unsigned int  							PowerCircuitCompMatrix::totalNumMarkings = 0;
PowerCircuitCompMatrix::NodeVector 		PowerCircuitCompMatrix::nodes;
unsigned int 							PowerCircuitCompMatrix::totalNumNodes = 0;
int										PowerCircuitCompMatrix::firstDeletedMarking = -1;
int										PowerCircuitCompMatrix::firstDeletedNode = -1;


PowerCircuitCompMatrix::PowerCircuitCompMatrix(int numInitNodes, int numInitCols)
:	matrix(numInitNodes,numInitCols),
	numTreedCols(0),
	numTreedNodes(0)
{
}


PowerCircuitCompMatrix::~PowerCircuitCompMatrix()
{
	for(int i = matrix.getNumRows() - 1; i >= 0; i--)
	{
		deleteNode(matrix.getRow(i));
	}


	// mark markings as deleted
	for(int j = matrix.getNumCols() - 1; j >= 0; j--)
	{
		for(std::list<Marking>::iterator markIter = matrix.colHdr(j).markings.begin(); markIter != matrix.colHdr(j).markings.end(); markIter++)
		{
			markings[markIter->id].deleted = true;
			markings[markIter->id].nextDeleted = firstDeletedMarking;
			firstDeletedMarking = markIter->id;
			totalNumMarkings--;
		}
	}

}


//-------------------------------------------------------------------------------------------
// Functions for checking integrity of nodes and markings

/**
 * checks if a marking is valid
 */
bool PowerCircuitCompMatrix::checkMarkingValid(const Marking& m) const
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
bool PowerCircuitCompMatrix::checkNodeValid(Node n) const
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
 * If a treed node is deleted, the user has to assure, that this does not change the property of being treed
 * for the other treed nodes.
 * Note that the values of other nodes and markings might be changed due to deletion of this node.
 */
void PowerCircuitCompMatrix::deleteNode(Row node)
{
	if(matrix.getRowNumber(node) < numTreedNodes) // matrix.rowHdr(node).treed
	{
		numTreedNodes--;
		if(matrix.getRowNumber(node) > 0)
			matrix.rowHdr(matrix.getRowNumber(node) - 1).BV = false;
	}

	Marking& mark = matrix.rowHdr(node).lambdaMarking;
	markings[mark.id].type = NORMAL; // necessary because successor markings cannot be deleted

	deleteMarking(mark);

	nodes[matrix.rowHdr(node).indexInNodes].deleted = true;
	nodes[matrix.rowHdr(node).indexInNodes].nextDeleted = firstDeletedNode;
	firstDeletedNode = matrix.rowHdr(node).indexInNodes;
	totalNumNodes--;

	matrix.deleteRow(node);
}


/**
 * deletes a columns of the matrix
 */
void PowerCircuitCompMatrix::deleteCol(Col col)
{
	assert(matrix.colHdr(col).markings.empty());
	if(matrix.getColNumber(col)<numTreedCols)
		numTreedCols--;
	matrix.deleteCol(col);
}


/**
 * Deletes the marking mark. If it is the last marking belonging to its column, the column is also deleted.
 * Only normal (not sucessor-) markings can be deleted.
 * Does NOT delete nodes which are not needed anymore.
 */
void PowerCircuitCompMatrix::deleteMarking(Marking& mark)
{
	assert(checkMarkingValid(mark));
	assert(markings[mark.id].type == NORMAL);
	Col col = markings[mark.id].col;

	// Find the marking in the list of its column and remove it from there.
	// True?: This is actually less efficient than the theoretical bound. True.
	for(std::list<Marking>::iterator markIter = matrix.colHdr(col).markings.begin(); markIter != matrix.colHdr(col).markings.end(); markIter++)
	{
		if(markIter->id == mark.id)
		{
			matrix.colHdr(col).markings.erase(markIter);
			break;
		}
	}

	if(matrix.colHdr(col).markings.empty())
		deleteCol(col);

	markings[mark.id].deleted = true;
	markings[mark.id].nextDeleted = firstDeletedMarking;
	totalNumMarkings--;
	firstDeletedMarking = mark.id;
	mark.id = -1;
	mark.pc = NULL;
}

void PowerCircuitCompMatrix::incMarkingRefCount(const Marking& mark)
{
	assert(checkMarkingValid(mark));
}
void PowerCircuitCompMatrix::decMarkingRefCount(Marking& mark)
{
	;
}

//-------------------------------------------------------------------------------------------
// Internal functions to create nodes and markings

/**
 * Creates a new marking pointing to col, updates markings and colHdr.
 * Does NOT create a new column.
 * Markings internally MUST NOT be created through other methods than this!!!
 */
Marking PowerCircuitCompMatrix::allocateNewMarking(Col col)
{
	int n;
	if(firstDeletedMarking == -1)
	{
		markings.resize(totalNumMarkings + 1);
		n = totalNumMarkings;
	}
	else
	{
		n = firstDeletedMarking;
		firstDeletedMarking = markings[n].nextDeleted;
	}

	markings[n].col = col;
	markings[n].type = NORMAL;
	markings[n].deleted = false;
	markings[n].pc = this;

	Marking mark =  Marking(this, n);
	matrix.colHdr(col).markings.insert(matrix.colHdr(col).markings.end(), mark);

	totalNumMarkings++;
	return mark;
}


/**
 * Takes a marking and turns it into the successor-marking of a new node.
 * The marking cannot be deleted anymore without deleting the node.
 * The newly created node is returned.
 * Nodes internally SHOULD NOT be created through other methods than this (except for cloning PCs)!!!
 */
Row PowerCircuitCompMatrix::newNodeFromMarking(const Marking& mark)
{
	assert(!markings[mark.id].deleted);
	int nodeNum;
	// insert node into static list of all nodes
	if(firstDeletedNode == -1)
	{
		nodes.resize(totalNumNodes + 1);
		nodeNum = totalNumNodes;
	}
	else
	{
		nodeNum = firstDeletedNode;
		firstDeletedNode = nodes[nodeNum].nextDeleted;
	}

	int n = matrix.getNumRows();
	matrix.insertZeroRow(n, RowHdr(mark, false, nodeNum));

	// make marking lambda
	markings[mark.id].type = LAMBDA;
	markings[mark.id].node = matrix.getRow(n);

	nodes[nodeNum].deleted = false;
	nodes[nodeNum].row = matrix.getRow(n);
	nodes[nodeNum].pc = this;
	totalNumNodes++;

	return matrix.getRow(n);
}


/**
 * Creates a new node with value one. If there is no treed node with value one, the new node is treed.
 */
Row PowerCircuitCompMatrix::newOneNode()
{
	Row node = newNodeFromMarking(newZeroMarking());
	if(numTreedNodes == 0)
	{
		insertCompactMarkingIntoTreed(matrix.rowHdr(node).lambdaMarking);
		moveNodeIntoTreedPartOfMatrix(node, 0);
	}
	return node;
}


/**
 * Clones a node, that means a new node is inserted with a copy of the successor-marking.
 * This successor-marking will be the only one pointing on its column in the matrix.
 * Even if the old node was treed, the new one is not!
 */
Row PowerCircuitCompMatrix::cloneNode(Row node)
{
	return newNodeFromMarking(newCopyColMarking(matrix.rowHdr(node).lambdaMarking));
}


/**
 * Creates a new marking with value one. I necessary a new node with value one is created.
 */
Marking PowerCircuitCompMatrix::newOneMarking()
{
	if(numTreedNodes == 0)
		newOneNode();
	int n = matrix.getNumCols();
	matrix.insertUnitCol(n, 0);
	return allocateNewMarking(matrix.getCol(n));
}


/**
 * Creates a new marking with value zero.
 */
Marking PowerCircuitCompMatrix::newZeroMarking()
{
	int n = matrix.getNumCols();
	matrix.insertZeroCol(n);
	return allocateNewMarking(matrix.getCol(n));
}


/**
 * Creates a new marking with the same value as mark. The column in the matrix is duplicated.
 */
Marking PowerCircuitCompMatrix::newCopyColMarking(const Marking& mark)
{
	assert(!markings[mark.id].deleted);
	int n = matrix.getNumCols();
	matrix.insertColCopy(markings[mark.id].col, n);
	return allocateNewMarking(matrix.getCol(n));
}


/**
 * Creates a new marking consisting of zeros and only one one at position onePos.
 */
Marking PowerCircuitCompMatrix::newUnitMarking(int onePos)
{
	int n = matrix.getNumCols();
	matrix.insertUnitCol(n,onePos);
	return allocateNewMarking(matrix.getCol(n));
}


//-------------------------------------------------------------------------------------------
// Auxiliary procedures to manipulate the way of storing markings in the matrix

/**
 * Ensures that the given marking is the only marking pointing on its column in the matrix.
 */
void PowerCircuitCompMatrix::separateMarkingFromCol(const Marking& m)
{
	Col col=markings[m.id].col;
	int n=matrix.getNumCols();
	if(matrix.colHdr(col).markings.size()>1)
	{
		assert(matrix.colHdr(n).markings.empty());
		matrix.colHdr(n).markings.push_back(m);

		matrix.insertColCopy(col,n);
		for(std::list<Marking>::iterator markIter=matrix.colHdr(col).markings.begin();markIter!=matrix.colHdr(col).markings.end();markIter++)
		{
			if(markIter->id==m.id)
			{
				matrix.colHdr(col).markings.erase(markIter);
				break;
			}
		}
		markings[m.id].col=matrix.getCol(n);
	}
}


/**
 * Changes markings pointing to col2 to point on targetCol, and then deletes the obsolete col2.
 */
void PowerCircuitCompMatrix::mergeCols(Col targetCol, Col col2)
{
	ColHdr& tColH=matrix.colHdr(targetCol);
	ColHdr& colH2=matrix.colHdr(col2);
	for(std::list<Marking>::iterator markIter=colH2.markings.begin(); markIter!=colH2.markings.end();markIter++)
	{
		markings[markIter->id].col=targetCol;
	}
	tColH.markings.insert(tColH.markings.end(), colH2.markings.begin(), colH2.markings.end());
	colH2.markings.clear();
	deleteCol(col2);
}



//-------------------------------------------------------------------------------------------
// Other auxiliary procedures dealing with treed nodes/markings

/**
 * Gives the Position of a treed node whose Lambda-marking is pointing on the given column col.
 * If there is no such node, the return value is -1.
 */
int PowerCircuitCompMatrix::getTreedNodePosToGivenCol(Col col)
{
	ColHdr colHdr=matrix.colHdr(col);
	for(std::list<Marking>::iterator markIter=colHdr.markings.begin() ; markIter!=colHdr.markings.end() ; markIter++)
	{
		if(markings[markIter->id].type==LAMBDA)
			if(matrix.getRowNumber(markings[markIter->id].node)<numTreedNodes)
				 return matrix.getRowNumber(markings[markIter->id].node);
	}
	return -1;
}


/**
 * Moves a node to the treed part of the matrix and marks it as treed.
 * NOTICE: Its Lambda-marking has to be treed already! Does NOT update the BV-vector!
 */
void PowerCircuitCompMatrix::moveNodeIntoTreedPartOfMatrix(Row node, int newPos)
{
	matrix.moveRow(node, newPos);
	numTreedNodes++;
}


/**
 * Returns the position, where a given compact column should be inserted in the treed part of the matrix.
 * If there is already a compact column with the same value, equal is set to true, otherwise it is set to false.
 */
int PowerCircuitCompMatrix::findNewPosOfCompactCol(Col col, bool& equal)
{
	int pos=0;
	Sign markEntry;
	int last=numTreedCols-1;
	equal=false;
	for(int i=numTreedNodes-1;i>=0;i--)
	{
		markEntry=matrix(i,col);
		while(pos<=last&&compareSigns(matrix(i,pos),markEntry)==-1)
			pos++;

		while(pos<=last&&compareSigns(matrix(i,last),markEntry)==1)
			last--;

		if(pos>last)
			return pos;
	}
	equal=true;
	return pos;
}


//-------------------------------------------------------------------------------------------
// Auxiliary procedures for incMarking

/**
 * Takes as input a column which is already compact or a compact column increased by one (and thus maybe not compact anymore)
 * Makes the column compact again. Does NOT move the column from or to the treed part of the PC.
 */
void PowerCircuitCompMatrix::compactifyFromBottom(Col col)
{
	int i=0;
	bool ready=false;
	if(matrix(0,col)==ZERO)
	{
		i++;
		if(!matrix.rowHdr(0).BV||matrix(1,col)==ZERO)//two Zeros at the beginning => is already compact
			return;
	}
	while(!ready)
	{
		if(!matrix.rowHdr(i).BV)
		{
			ready=true;
		}
		else if(matrix(i+1,col)==ZERO) // in this case node i is NOT the last node in tree-rep
		{
			ready=true;
		}
		else //not ready
		{
			if(matrix(i+1,col)==MINUS && matrix(i,col)==PLUS)
			{
				matrix.setEntry(i+1,col,ZERO);
				matrix.setEntry(i,col,MINUS);
				ready=true;
			}
			else  //in this case is (matrix(i+1,col)==PLUS && matrix(i,col)==PLUS) -- not ready
			{
				assert(matrix(i+1,col)==PLUS && matrix(i,col)==PLUS);

				if(!matrix.rowHdr(i+1).BV){
					insertNewPowerOfTwoNode(i+2);
				}
				matrix.setEntry(i+2,col,PLUS); //was before ==0 (due to compactness)
				matrix.setEntry(i+1,col,ZERO);
				matrix.setEntry(i,col,MINUS);
				i+=2;
			}
		}
	}
}


/**
 * Creates a new node with value 2^power (hence with successor-marking of value power) and inserts it into the
 * treed part of the PC. All nodes with values 2^i (0 <= i < power) have to exist already!!!
 */
void PowerCircuitCompMatrix::insertNewPowerOfTwoNode(unsigned int power)
{
	if(matrix.rowHdr(power-1).BV)
		return;
	Row node=newNodeFromMarking(newZeroMarking());
	Sign temp[32];

	unsigned int numReqNodes=calculateCompactRepresentation(power,temp);
	for(unsigned int i=0;i<numReqNodes;i++)
	{
		matrix.setEntry(i,markings[matrix.rowHdr(node).lambdaMarking.id].col,temp[i]);
	}

	insertCompactMarkingIntoTreed(matrix.rowHdr(node).lambdaMarking);
	moveNodeIntoTreedPartOfMatrix(node,power);
	matrix.rowHdr(power-1).BV=true;

	if(power<numTreedNodes-1)
	{
		numReqNodes=calculateCompactRepresentation(power+1,temp);
		bool equal=true;
		Col compareCol=markings[matrix.rowHdr(power+1).lambdaMarking.id].col;

		for(unsigned int i=0;i<numReqNodes;i++)
		{
			if(matrix(i,compareCol)!=temp[i])
			{
				equal=false;
				break;
			}
		}
		if(equal)
			for(unsigned int i=numReqNodes;i<power;i++)
			{
				if(matrix(i,compareCol)!=ZERO)
				{
					equal=false;
					break;
				}
			}
		matrix.rowHdr(power).BV=equal;
	}
	else
	{
		matrix.rowHdr(power).BV=false;
	}
}


/**
 * Calculates the compact Representation of a given integer.
 */
unsigned int PowerCircuitCompMatrix::calculateCompactRepresentation(int n, Sign result[32])
{
	unsigned int temp=(n>=0)?n:(-n);
	int uebertrag=0;
	unsigned int res=0;
	for(int i=0;i<32;i++)
	{
		if(uebertrag+(temp&1)==0)
		{
			result[i]=ZERO;
			uebertrag=0;
		}
		else if(uebertrag+(temp&1)==2)
		{
			result[i]=ZERO;
			uebertrag=1;
		}
		else if((temp&2)!=0)
		{
			result[i]=MINUS;
			uebertrag=1;
		}
		else
		{
			result[i]=PLUS;
			res=i+1;
			uebertrag=0;
		}
		temp=temp>>1;
	}
	if(n<0)
	{
		for(int i=0;i<32;i++)
		{
			if(result[i]==MINUS)
				result[i]=PLUS;
			else if(result[i]==PLUS)
				result[i]=MINUS;
		}
	}
	assert(0<=res&&res<=32);
	return res;
}


//------------------------------------------------------------------------------------------------
// Auxiliary procedures for reduce()

/**
 * Sets the BV-value of a node correctly. The node has to be in the treed part of the PC.
 */
void PowerCircuitCompMatrix::setBV(Row node)
{
	unsigned int nodePos=matrix.getRowNumber(node);
	if(nodePos<numTreedNodes-1)
	{
		Marking lambdaPlusOne=incMarking(matrix.rowHdr(node).lambdaMarking);
		matrix.rowHdr(node).BV = matrix.columnsEqual(markings[matrix.rowHdr(nodePos+1).lambdaMarking.id].col, markings[lambdaPlusOne.id].col);
		deleteMarking(lambdaPlusOne);
	}
	else
	{
		matrix.rowHdr(node).BV=false;
	}
}


/**
 * Creates a new node with the double value of a treed node "node" and inserts it into the treed part of the PC.
 */
void PowerCircuitCompMatrix::newDoubleNode(Row node)
{
	if(matrix.getRowNumber(node)>=numTreedNodes)
	{
		std::cout<<"Error in newDoubleNode: node not treed! RowIndex:"<<matrix.getRowNumber(node)<<" numTreedNodes: "<<numTreedNodes<<std::endl;
		assert(false);
		return;
	}
	if(matrix.rowHdr(node).BV) // if there is alreay a node with the double value
	{
		return;
	}

	Marking doubleMarking=incMarking(matrix.rowHdr(node).lambdaMarking);
	insertCompactColIntoTreed(markings[doubleMarking.id].col);
	Row newNode=newNodeFromMarking(doubleMarking);

	moveNodeIntoTreedPartOfMatrix(newNode, matrix.getRowNumber(node)+1);
	matrix.rowHdr(node).BV=true;
	setBV(newNode);
}


/**
 * Makes a column compact. Does Not move it to the treed part of the Circuit.
 */
void PowerCircuitCompMatrix::compactify(Col col)
{
	unsigned int i=0;
	Sign carry=ZERO;
	Sign tmpResult;
	Sign tmpCarry;
	while(i<numTreedNodes)
	{
		addSigns(matrix(i,col),carry,tmpResult,tmpCarry);
		carry=tmpCarry;
		matrix.setEntry(i,col,tmpResult);

		if(!matrix.rowHdr(i).BV||matrix(i,col)==ZERO)
		{
			i++;
		}
		else// i is NOT the last node in tree-rep
		{
			if(matrix(i+1,col)==ZERO)
			{
				i+=2;
			}
			else
			{
				if(matrix(i+1,col)==MINUS&&matrix(i,col)==PLUS)
				{
					matrix.setEntry(i+1,col,ZERO);
					matrix.setEntry(i,col,MINUS);
					i+=2;
				}
				else if(matrix(i+1,col)==PLUS&&matrix(i,col)==MINUS)
				{
					matrix.setEntry(i+1,col,ZERO);
					matrix.setEntry(i,col,PLUS);
					i+=2;
				}
				else if(!((matrix(i+1,col)==PLUS&&matrix(i,col)==PLUS)||((matrix(i+1,col)==MINUS&&matrix(i,col)==MINUS))))
				{
					std::cout<<(int)matrix(i+1,col)<<", "<<(int)matrix(i,col)<<std::endl;
					assert(false);
				}
				else
				{
					tmpResult=matrix(i+1,col);
					if(!matrix.rowHdr(i+1).BV){
						newDoubleNode(matrix.getRow(i+1));
						matrix.setEntry(i+2,col,tmpResult);
					}
					else //if(nodes[nodeOrder[i+1]].BV)
					{
						carry=tmpResult;
					}
					matrix.setEntry(i+1,col,ZERO);
					matrix.setEntry(i,col,negateSign(tmpResult));
					i+=2;
				}
			}
		}
	}
	if(carry!=ZERO)
	{
		newDoubleNode(matrix.getRow(numTreedNodes-1));
		matrix.setEntry(numTreedNodes-1,col,carry);
	}
}


/**
 * Checks all non-treed columns if they use oldNode, and, if so, oldNode is replaced by the node newNode.
 * newNode has to be in the treed part of the PC!!! If a marking uses both oldNode and newNode, new nodes are created
 * and inserted into the treed part of the PC
 */
void PowerCircuitCompMatrix::removeDoubleNodesFromMarkings(Row oldNode, Row newNode)
{
	Sign result;
	Sign carry;
	for(unsigned int j=numTreedCols;j<matrix.getNumCols();j++)
	{
		if(matrix(oldNode,j)!=ZERO)
		{
			addSigns(matrix(oldNode,j),matrix(newNode,j),result,carry);
			if(carry!=0)
			{
				unsigned int i= matrix.getRowNumber(newNode);
				while(carry!=0)
				{
					assert(i<numTreedNodes);
					if(matrix(i,j)==carry)
					{
						if(!matrix.rowHdr(i).BV)
							newDoubleNode(matrix.getRow(i));
						matrix.setEntry(i,j,ZERO);
						i++;
					}
					else
					{
						Sign tmpCarry=carry;
						addSigns(matrix(i,j),tmpCarry,result,carry);
						matrix.setEntry(i,j,result);
						assert(carry==ZERO);
					}
				}
			}
			else
			{
				matrix.setEntry(newNode,j, result);
			}
			matrix.setEntry(oldNode,j,ZERO);
		}
	}
}


/**
 * 	Inserts the column of the marking mark into the treed part of the PC. The column has to be compact already!!
 * 	Returns the new position of the column.
 */
int PowerCircuitCompMatrix::insertCompactMarkingIntoTreed(const Marking& mark)
{
	assert(checkMarkingValid(mark));
	return insertCompactColIntoTreed(markings[mark.id].col) ;
}


/**
 * 	Inserts the column col into the treed part of the PC. The column has to be compact already!!
 * 	Returns the new position of the column.
 */
int PowerCircuitCompMatrix::insertCompactColIntoTreed(Col col)
{
	assert(matrix.getColNumber(col)>=(unsigned int)numTreedCols);
	bool equal;
	unsigned int newPos=findNewPosOfCompactCol(col, equal);

	if(newPos==numTreedCols||!equal)
	{
		matrix.moveCol(col, newPos);
		numTreedCols++;
	}
	else
	{
		mergeCols(matrix.getCol(newPos),col);
	}
	return newPos;
}


/**
 * Inserts the node into the treed part of the PC.
 * Its successor-marking has to be in the treed part of the PC already!!
 */
int PowerCircuitCompMatrix::insertNodeIntoTreed(Row node)
{
	Col markingCol=markings[matrix.rowHdr(node).lambdaMarking.id].col;
	int newPos= getTreedNodePosToGivenCol(markingCol);

	if(newPos==-1)// if there is no treed node with the same Lambda-Marking as node
	{
		//determine newPos
		unsigned int i=0;
		while(i<numTreedNodes&&matrix.getColNumber(markings[matrix.rowHdr(i).lambdaMarking.id].col)<matrix.getColNumber(markingCol))
			i++;
		newPos=i;

		moveNodeIntoTreedPartOfMatrix(node, newPos);
		if(newPos>0)
			setBV(matrix.getRow(newPos-1));
		setBV(matrix.getRow(newPos));
	}
	else // if there is already a treed node with the same Lambda-Marking
	{
		removeDoubleNodesFromMarkings(node, matrix.getRow(newPos));
		deleteNode(node);
	}
	return newPos;
}


/**
 * Moves the columns in colList to the treed part of the PC.
 * These columns have to have their support in the treed part of the PC!!
 * colList is emptied in this process
 */
void PowerCircuitCompMatrix::moveColsToTreed(std::list<Col>& colList)
{
	while(!colList.empty())
	{
		compactify(colList.front());
		insertCompactColIntoTreed(colList.front());
		colList.pop_front();
	}
}


/**
 * Recursive procedure to sort topologically starting with node nodesToSort[nodeIndex].
 * This procedure has to be called for all nodes which should be sorted topologically.
 */
void PowerCircuitCompMatrix::topSortNode(unsigned int nodeIndex, std::vector<Row>& nodesToSort, std::vector<bool>& visited, std::vector<Row>& topSortPerm)
{
	Row tmpNode=nodesToSort[nodeIndex];
	visited[nodeIndex]=true;

	Col col=markings[matrix.rowHdr(tmpNode).lambdaMarking.id].col;
	for(unsigned int i = 0 ; i<nodesToSort.size() ; i++)
	{
		if((!visited[i])&&matrix(nodesToSort[i], col)!=ZERO)
			topSortNode(i, nodesToSort, visited, topSortPerm);
	}

	topSortPerm.push_back(tmpNode);
}


/**
 *  Moves all nodes in nodeList and columns in colList to the treed part of the PC.
 *  At the end nodeList is the empty list!!!
 *  All columns of successor-markings of nodes in nodeList have to be in colList.
 *  The support of columns in colList has to be in nodeList or the already treed part of the PC.
 */
void PowerCircuitCompMatrix::extendTree(std::vector<Row>& nodeList, std::vector<Col>& colList)
{
	std::vector<Row>  topSortPerm;
	int numNodesToSort=nodeList.size();
	topSortPerm.reserve(numNodesToSort);
	std::vector<bool>  visited;
	visited.assign(numNodesToSort, false);

	//sort topologically
	for(int i = 0 ; i<numNodesToSort ; i++)
	{
		if((!visited[i]))
			topSortNode(i, nodeList, visited, topSortPerm);
	}

	assert((int)topSortPerm.size()==numNodesToSort);
	std::list<Col> treedSuppCols;
	std::vector<std::list<Col> > lastNodeInColList(numNodesToSort);

	//sort columns by the last node in topSortPerm in their support
	bool treed;
	for(std::vector<Col>::iterator colIter=colList.begin() ; colIter!=colList.end() ; colIter++ )
	{
		treed=true;
		for(int i=numNodesToSort-1; i>=0 ; i--)
		{
			if(matrix(topSortPerm[i],*colIter)!=ZERO)
			{
				lastNodeInColList[i].push_back(*colIter);
				treed=false;
				break;
			}
		}
		if(treed)
			treedSuppCols.push_back(*colIter);
	}

	moveColsToTreed(treedSuppCols);

	for(int i=0; i < numNodesToSort ; i++)
	{
		insertNodeIntoTreed(topSortPerm[i]);
		moveColsToTreed(lastNodeInColList[i]);
	}
}


//-------------------------------------------------------------------------------------------------
// implementation of private abstract methods of PowerCircuit

/**
 * Creates a new Marking with the inverse of mark.
 */
Marking PowerCircuitCompMatrix::invMarking(const Marking& mark)
{
	assert(checkMarkingValid(mark));
	int n=matrix.getNumCols();
	matrix.insertInvCol(markings[mark.id].col,n);
	return allocateNewMarking(matrix.getCol(n));
}


/**
 * Creates a new Marking with the sum of mark1 and mark2. The user has to assure that for each node,
 * the sum lies in the range [-1,+1]!!!
 */
Marking PowerCircuitCompMatrix::addMarkings(const Marking& mark1, const Marking& mark2)
{
	assert(checkMarkingValid(mark1));
	assert(checkMarkingValid(mark2));
	int n=matrix.getNumCols();
	matrix.insertColSum(markings[mark1.id].col, markings[mark2.id].col,n);
	return allocateNewMarking(matrix.getCol(n));
}


/**
 * Creates a new Marking m with the intersection of mark1 and mark2, that is
 * m(p) = mark1(p), 	if (mark1(p)==mark2(p))
 * m(p) = 0 			else
 * for every node p.
 */
Marking PowerCircuitCompMatrix::intersectMarkings(const Marking& mark1, const Marking& mark2)
{
	assert(checkMarkingValid(mark1));
	assert(checkMarkingValid(mark2));
	int n=matrix.getNumCols();
	matrix.insertIntersectionCol(markings[mark1.id].col, markings[mark2.id].col,n);
	return allocateNewMarking(matrix.getCol(n));
}


/**
 * Increases given treed marking by one. The result is a compact but NOT treed marking
 */
Marking PowerCircuitCompMatrix::incMarking(const Marking& mark)
{
	assert(checkMarkingValid(mark));

	if(matrix.getColNumber(markings[mark.id].col)>=numTreedCols) //if the marking is not reduced
	{
		Marking result=newCopyColMarking(mark);
		if(numTreedNodes==0||matrix(0,markings[mark.id].col)==PLUS)
		{
			Row oneNode = newOneNode();
			matrix.setEntry(oneNode,markings[result.id].col,PLUS);
		}
		else if(matrix(0,markings[mark.id].col)==ZERO)
		{
			matrix.setEntry(0,markings[result.id].col,PLUS);
		}
		else //(matrix(0,col)==MINUS)
		{
			matrix.setEntry(0,markings[result.id].col,ZERO);
		}
		return result;
	}
	else		//if the marking is reduced
	{
		Marking result=newCopyColMarking(mark);
		if(numTreedNodes==0)
			newOneNode();
		Row oneNode= matrix.getRow(0);

		//Increment the value of the result column by one
		Sign matrixEntry = matrix(oneNode,markings[result.id].col);

		if(matrixEntry==MINUS)
		{
			matrix.setEntry(oneNode,markings[result.id].col,ZERO);
		}
		else if(matrixEntry==ZERO)
		{
			matrix.setEntry(oneNode,markings[result.id].col,PLUS);
		}
		else //matrixEntry==PLUS
		{
			Row twoNode;
			if(matrix.rowHdr(oneNode).BV==true)
			{
				twoNode=matrix.getRow(1);
			}
			else 	//can only occur if the tree-representation has only one node
			{
				assert(numTreedNodes==1);
				twoNode=newNodeFromMarking(newUnitMarking(matrix.getRowNumber(oneNode)));

				insertCompactMarkingIntoTreed(matrix.rowHdr(twoNode).lambdaMarking);
				moveNodeIntoTreedPartOfMatrix(twoNode, 1);
				matrix.rowHdr(oneNode).BV=true;
			}
			matrix.setEntry(oneNode,markings[result.id].col,ZERO);
			matrix.setEntry(twoNode,markings[result.id].col,PLUS);
		}

		compactifyFromBottom(markings[result.id].col);
		return result;
	}
}


/**
 * Clones a marking. That means all nodes in the support of the marking are cloned and a new marking with
 * the same value as the original one, using only the newly created nodes, is created
 */
Marking PowerCircuitCompMatrix::cloneMarking(const Marking& mark)
{
	assert(checkMarkingValid(mark));
	Marking result=newZeroMarking();

	int oldNumNodes=matrix.getNumRows();
	Col clonedCol=markings[result.id].col;

	for(int i=0;i<oldNumNodes;i++)
	{
		if(matrix(i,markings[mark.id].col)!=ZERO){
			Row node=cloneNode(matrix.getRow(i));
			matrix.setEntry(node,clonedCol,matrix(i,markings[mark.id].col));
		}
	}
	return result;
}


/**
 * Clones a node, that means a new node is inserted with a copy of the successor-marking.
 * This successor-marking will be the only one pointing on its column in the matrix.
 * Even if the old node was treed, the new one is not!
 */
Node PowerCircuitCompMatrix::cloneNode(Node n)
{
	assert(checkNodeValid(n));
	Row newNode=cloneNode(nodes[n.id].row);
	return Node(this,matrix.rowHdr(newNode).indexInNodes);
}


/*
 * Return a deep copy of a marking (that means a new marking with the same nodes as the old one).
 */
Marking	PowerCircuitCompMatrix::copyMarking(const Marking& m)
{
	assert(checkMarkingValid(m));
	return newCopyColMarking(m);
}


/**
 * Return whether a marking is in the reduces part of the PC.
 */
bool PowerCircuitCompMatrix::isMarkingReduced(const Marking& m) const
{
	assert(checkMarkingValid(m));
	return (matrix.getColNumber(markings[m.id].col)<numTreedCols);
}


/**
 * Compares the values of two markings. Returns
 * -1, if m1<m2,
 * 0,  if m1==m2
 * 1,  if m1>m1
 */
int	PowerCircuitCompMatrix::compareMarkings(const Marking& m1, const Marking& m2)
{
	assert(checkMarkingValid(m1));
	assert(checkMarkingValid(m2));
	if(matrix.getColNumber(markings[m1.id].col)>=numTreedCols||matrix.getColNumber(markings[m2.id].col)>=numTreedCols)
	{
		throw "Only markings in the reduced part of the circuit can be compared.";
	}
	Col col1=markings[m1.id].col;
	Col col2=markings[m2.id].col;
	for(int i=numTreedNodes-1;i>=0;i--)
	{
		int d=compareSigns(matrix(i,col1),matrix(i,col2));
		if(d!=0)
			return d;
	}
	return 0;
}


/**
 *	Returns the index of a node in the reduces part of the PC
 *	(thus the number of smaller nodes in the reduced part).
 *	If the node is not in the reduced part, -1 is returned.
 */
int	PowerCircuitCompMatrix::getRedNodeOrd(Node n)
{
	assert(checkNodeValid(n));
	if(matrix.getRowNumber(nodes[n.id].row)<numTreedNodes)
		return matrix.getRowNumber(nodes[n.id].row);
	else
		return -1;
}


/**
 * Returns the Sign of a given node n in a marking m. 1 for PLUS, -1 for MINUS.
 */
int	PowerCircuitCompMatrix::getNodeSignInMarking(Node n, const Marking& m) const
{
	assert(checkNodeValid(n));
	assert(checkMarkingValid(m));
	Sign sign=matrix(nodes[n.id].row,markings[m.id].col);
	if(sign==ZERO)
		return 0;
	else if(sign==MINUS)
		return -1;
	else
		return 1;
}


/**
 * Returns whether the marking m is a successor-marking.
 */
bool PowerCircuitCompMatrix::isSuccessorMarking(const Marking& m) const
{
	assert(checkMarkingValid(m));
	return (markings[m.id].type==LAMBDA);
}


/**
 * Takes as input a marking m. If m is a successor-marking,
 * the node whose successor-marking it is is returned.
 * Otherwise the undefined node is returned.
 */
Node PowerCircuitCompMatrix::getIncidentNode(const Marking& m)
{
	assert(checkMarkingValid(m));
	if(markings[m.id].type==LAMBDA)
		return Node(this,matrix.rowHdr(markings[m.id].node).indexInNodes);
	else
		return Node(this,-1);
}


/**
 * Returns the smallest non-zero node of a treed marking m.
 * If m is the zero-marking the undefined node is returned.
 */
Node PowerCircuitCompMatrix::getSmallestNodeInMarking(const Marking& m)
{
	assert(checkMarkingValid(m));
	Col col=markings[m.id].col;
	if(matrix.getColNumber(col)>=numTreedCols)
		throw "getSmallestNodeInMarking: Marking has to be treed!";
	for(unsigned int i=0;i<numTreedNodes;i++)
	{
		if(matrix(i,col)!=ZERO)
			return Node(this,matrix.rowHdr(i).indexInNodes);
	}
	return Node(this,-1);
}


/**
 * Returns the successor-marking of the node n.
 */
Marking	PowerCircuitCompMatrix::getSuccMarking(Node n)
{
	assert(checkNodeValid(n));
	return matrix.rowHdr(nodes[n.id].row).lambdaMarking;
}



//---------------------------------------------------------------------------------------------
// public methods

/**
 * Reduces the whole PC and all its markings
 */
void PowerCircuitCompMatrix::reduce()
{
	std::vector<Row>  nodeList;
	nodeList.reserve(matrix.getNumRows()-numTreedNodes);
	std::vector<Col> colList;
	colList.reserve(matrix.getNumCols()-numTreedCols);

	//put nodes, which should be treed (i.e. all non-treed nodes) into temporary list
	for(unsigned int i=numTreedNodes;i<matrix.getNumRows();i++)
	{
		nodeList.push_back(matrix.getRow(i));
	}
	for(unsigned int j = numTreedCols ; j<matrix.getNumCols() ; j++)
	{
		colList.push_back(matrix.getCol(j));
	}

	extendTree(nodeList, colList);
}


void PowerCircuitCompMatrix::reduce(std::vector<Node>& nodeVector, std::vector<Marking>& markingVector)
{
	std::vector<Row> nodeRowList;
	nodeRowList.reserve(nodeVector.size());
	std::vector<Col> colList;
	colList.reserve(matrix.getNumCols()-numTreedCols);
	std::vector<bool> isInNodeList(matrix.getNumRows()-numTreedNodes);
	for(int i = 0; i<(int)matrix.getNumRows()-(int)numTreedNodes; i++)
	{
		isInNodeList[i]=false;
	}

	for(std::vector<Node>::iterator nodeIter=nodeVector.begin() ; nodeIter!=nodeVector.end(); nodeIter++)
	{
		if(matrix.getRowNumber(nodes[nodeIter->id].row)>=numTreedNodes)
			isInNodeList[matrix.getRowNumber(nodes[nodeIter->id].row)-numTreedNodes]=true;
	}


	//prepare nodeRowList (node to be treed) and exclude-list (nodes which should not be treed) for extendTree
	std::list<Row> excludeList;
	std::vector<bool> colNeeded;
	colNeeded.assign(matrix.getNumCols()-numTreedCols, false);

	for(int i = 0; i<(int)matrix.getNumRows()-(int)numTreedNodes; i++)
	{
		if(!isInNodeList[i])
		{
			excludeList.push_back(matrix.getRow(i+numTreedNodes));
		}
		else
		{
			nodeRowList.push_back(matrix.getRow(i+numTreedNodes));
			colNeeded[matrix.getColNumber(markings[matrix.rowHdr(i+numTreedNodes).lambdaMarking.id].col)+numTreedCols] = true;
		}
	}
#ifdef DEBUG
	//////////////////////////////////////////////
	//security check
	for(std::vector<Row>::iterator nodeIter=nodeRowList.begin() ; nodeIter!=nodeRowList.end(); nodeIter++)
	{
		Col col=markings[matrix.rowHdr(*nodeIter).lambdaMarking.id].col;
		for(std::list<Row>::iterator nodeIter = excludeList.begin(); nodeIter!=excludeList.end() ; nodeIter++)
		{
			if(matrix(*nodeIter,col)!=ZERO)
			{
				throw "Invalid nodeList for reduce()!!";
			}
		}
	}
	/////////////////////////////////////////////
#endif

	for(std::vector<Marking>::iterator markIter=markingVector.begin() ; markIter!=markingVector.end(); markIter++)
	{
		colNeeded[matrix.getColNumber(markings[markIter->id].col)+numTreedCols] = true;
	}

	for(unsigned int j = numTreedCols ; j<matrix.getNumCols() ; j++)
	{
		if(colNeeded[j - numTreedCols])
			colList.push_back(matrix.getCol(j));
	}

	extendTree(nodeRowList, colList);
}


/**
 * Creates a new marking with value val. The therefore required nodes are inserted into the treed part of the PC.
 */
Marking	PowerCircuitCompMatrix::createMarking(int val)
{
	if(val==0)
		return newZeroMarking();
	else if(val==1)
		return newOneMarking();
	else
	{
		Marking m =newZeroMarking();
		Col col= markings[m.id].col;

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
		if(numTreedNodes==0)
			newOneNode();
		if((temp&0x1)==1)
			matrix.setEntry(i,col,sign);
		temp=temp>>1;
		i++;

		while(temp!=0)
		{

			assert((unsigned int)i<8*sizeof(int));
			if(!matrix.rowHdr(i-1).BV)
				insertNewPowerOfTwoNode(i);

			if((temp&0x1)==1)
				matrix.setEntry(i,col,sign);
			temp=temp>>1;
			i++;
		}
		return m;
	}
}


/**
 *  Creates a new marking, with the nodes in nodeList set to 1.
 */
Marking	PowerCircuitCompMatrix::createMarking(const std::list<Node>& nodeList)
{
	Marking mark=newZeroMarking();
	Col col=markings[mark.id].col;
	for(std::list<Node>::const_iterator nodeIter = nodeList.begin() ; nodeIter != nodeList.end() ; nodeIter++)
	{
		matrix.setEntry(nodes[nodeIter->id].row,col,PLUS);
	}
	return mark;
}


/**
 *  Creates a new marking, with the nodes in nodeList set to 1.
 */
Marking	PowerCircuitCompMatrix::createMarking(const std::list<Node>& nodeList, const std::list<int>& signList)
{
	Marking mark=newZeroMarking();
	Col col=markings[mark.id].col;
	std::list<int>::const_iterator signIter = signList.begin();
	for(std::list<Node>::const_iterator nodeIter = nodeList.begin() ;	nodeIter != nodeList.end() &&  signIter != signList.end(); nodeIter++)
	{
		if(*signIter > 0)
			matrix.setEntry(nodes[nodeIter->id].row,col,PLUS);
		else
			matrix.setEntry(nodes[nodeIter->id].row,col,MINUS);

		signIter++;
	}
	return mark;
}


/**
 * Creates a new node with successor-marking m. The new successor-marking points on the same column as m.
 */
Node PowerCircuitCompMatrix::createNode(const Marking& m)
{
	assert(checkMarkingValid(m));
	Marking mark =allocateNewMarking(markings[m.id].col);
	Row node=newNodeFromMarking(mark);
	return Node(this,matrix.rowHdr(node).indexInNodes);
}


void PowerCircuitCompMatrix::checkCyclesRecursive(Row node,std::vector<bool>  visited)
{
	if(visited[matrix.getRowNumber(node)])
	{
		draw("Cyclegraph");
		throw "Cycle in Graph";
	}
	visited[matrix.getRowNumber(node)] = true;
	Col col = markings[matrix.rowHdr(node).lambdaMarking.id].col;
	for(unsigned int i = 0 ; i < matrix.getNumRows() ; i++)
	{
		if(matrix(i,col) != ZERO)
			checkCyclesRecursive(matrix.getRow(i), visited);
	}
	visited[matrix.getRowNumber(node)] = false;
}

void PowerCircuitCompMatrix::checkCycles()
{
	std::vector<bool> visited;
	visited.assign(matrix.getNumRows(),false);
	for(unsigned int i = 0 ; i < matrix.getNumRows() ; i++)
	{
		checkCyclesRecursive(matrix.getRow(i), visited);
	}
}


/**
 * Adds the marking  p to all successor-markings of nodes in the support of \e m. (So, if supp(m) and supp(p) are
 * disjoint, that means, that arrows are drawn from all nodes in supp(m) to p).
 * The user has to assure that, for each node, the sum lies in the range [-1,+1] (for example using clone())!!!
 * The user is required to check that the resulting graph has no cycles.
 * All nodes in supp(m) have to be non-treed. (use clone() to assure that).
 */
void PowerCircuitCompMatrix::connect(const Marking& m, const Marking& p)
{
	assert(checkMarkingValid(m));
	assert(checkMarkingValid(p));
	Col tcol=markings[m.id].col;

	for(unsigned int i=0;i<matrix.getNumRows();i++)
	{
		if(matrix(i,markings[m.id].col)!=ZERO)
		{
			if(i<numTreedNodes)
			{
				throw "Can only connect non-treed nodes.";
			}
			separateMarkingFromCol(matrix.rowHdr(i).lambdaMarking);
			matrix.addToCol(markings[p.id].col,markings[matrix.rowHdr(i).lambdaMarking.id].col);
		}
	}

#ifdef DEBUG
 	checkCycles();
#endif
}

/**
 * The same as connect, only that the marking q is subtracted instead of added.
 */
void PowerCircuitCompMatrix::connectInv(const Marking& m, const Marking& q)
 {
 	assert(checkMarkingValid(m));
 	assert(checkMarkingValid(q));
 	Col tcol=markings[m.id].col;
 	Marking p=invMarking(q);

 	for(unsigned int i=0;i<matrix.getNumRows();i++)
 	{
 		if(matrix(i,markings[m.id].col)!=ZERO)
 		{
 			if(i<numTreedNodes)
 			{
 				throw "Can only connect non-treed nodes.";
 			}
 			separateMarkingFromCol(matrix.rowHdr(i).lambdaMarking);
 			matrix.addToCol(markings[p.id].col,markings[matrix.rowHdr(i).lambdaMarking.id].col);
 		}
 	}
 	deleteMarking(p);
#ifdef DEBUG
 	checkCycles();
#endif
 }


/**
 * Returns the node with index ord (that means there are ord smaller nodes in the reduced part of the PC)
 * If there are less than ord node, the undefined node is returned
 */
Node PowerCircuitCompMatrix::getReducedNode(unsigned int ord)
{
	if(ord<numTreedNodes)
		return Node(this,matrix.rowHdr(ord).indexInNodes);
	else
		return Node(this,-1);
}


/**
 * Returns a list of all nodes in the PowerCircuit.
 */
std::list<Node>	PowerCircuitCompMatrix::getNodes()
{
	std::list<Node> nodes;
	for(unsigned int i=0;i<matrix.getNumRows();i++)
	{
		nodes.push_back(Node(this,matrix.rowHdr(i).indexInNodes));
	}
	return nodes;
}


/**
 * Returns a list of all markings in the PowerCircuit.
 */
std::list<Marking>	PowerCircuitCompMatrix::getMarkings()
{
	std::list<Marking> marks;
	for(unsigned int j=0;j<matrix.getNumCols();j++)
	{
		marks.insert(marks.end(),matrix.colHdr(j).markings.begin(),matrix.colHdr(j).markings.end());
	}
	return marks;
}


/**
 * Returns a list of all nodes which are non-zero in the marking m.
 */
std::list<Node>	PowerCircuitCompMatrix::getMarkingNodes(const Marking& m)
{
	assert(checkMarkingValid(m));
	std::list<Node> nodes;
	Col col=markings[m.id].col;
	for(unsigned int i=0;i<matrix.getNumRows();i++)
	{
		if(matrix(i,col)!=ZERO)
			nodes.push_back(Node(this,matrix.rowHdr(i).indexInNodes));
	}
	return nodes;
}


/**
 * Prints the matrix of the PC.
 */
void PowerCircuitCompMatrix::print(std::ostream& os)
{
	os<<matrix.getNumRows()<<" x "<<matrix.getNumCols()<<" Matrix:"<<std::endl;

	os<<"totalNumMarkings: "<<totalNumMarkings<<std::endl;
	os << std::setw(4) << " ";
	for(unsigned int j = 0; j < matrix.getNumCols(); j++){
		if(j<numTreedCols)
			os << std::setw(4) << "t";
		else
			os << std::setw(4) << " ";
	}
	os << std::endl;

	for(unsigned int i = 0; i < matrix.getNumRows(); i++)
	{
		if(i<numTreedNodes)
			os << std::setw(3) << "t"<< (matrix.rowHdr(i).BV?1:0);
		else
			os << std::setw(4) << " ";
		for(unsigned int j = 0; j <matrix.getNumCols(); j++){
			os << std::setw(4) << signToString(matrix(i,j));
		}
		os << std::endl;
	}
	unsigned int maxNumMarkingsInCol=0;
	for(unsigned int j = 0; j < matrix.getNumCols(); j++)
	{
		 maxNumMarkingsInCol=std::max( maxNumMarkingsInCol,(unsigned int)matrix.colHdr(j).markings.size());
	}
	for(unsigned int k = 0; k <maxNumMarkingsInCol; k++)
	{
		os << std::setw(4) << " ";
		for(unsigned int j = 0; j < matrix.getNumCols(); j++){
			if(matrix.colHdr(j).markings.size()>k)
			{
				os.fill('0');

				//find the k-th element in the list
				std::list<Marking>::iterator markIter=matrix.colHdr(j).markings.begin();
				for(unsigned int s = 0 ; s < k ; s++)
					markIter++;

				if(markings[markIter->id].type==LAMBDA)
					os << " L"<< std::setw(2) <<matrix.getRowNumber(markings[markIter->id].node);
				else
					os << " M"<< std::setw(2) <<markIter->id;
				os.fill(' ');
			}
			else
				os << std::setw(4) << " ";
		}
		os<<std::endl;
	}

	os << std::endl;
	os.flush();
}


/**
 * Deletes all nodes of the marking. (So if other markings use those
 * nodes, their values may be changed!!!)
 */
 void PowerCircuitCompMatrix::remove(const Marking& m)
 {
	 assert(checkMarkingValid(m));
	 std::list<Node> nodeList=getMarkingNodes(m);
	 while(!nodeList.empty())
	 {
		 Node node = nodeList.front();
		 nodeList.pop_front();
		 deleteNode(nodes[node.id].row);
	 }
 }


 /**
  * Makes a copy of a PowerCircuit. Only the markings in the markingsToClone list are copied.
  * Treed nodes and markings are kept treed.
  * Has to be deallocated via delete, once it is not needed anymore!!!
  */
 PowerCircuit* PowerCircuitCompMatrix::clone(std::vector<Marking>& markingsToKeep)
 {
	PowerCircuitCompMatrix* result = new PowerCircuitCompMatrix();
	int numCols=matrix.getNumCols();
	int numRows=matrix.getNumRows();

	//write in a vector which columns should be copied
	std::vector<bool> keepCol(numCols);
	for(int j = 0 ; j < numCols ; j++)
		keepCol[j]=false;
	for(int i=0 ; i < numRows ; i++)
	{
		keepCol[matrix.getColNumber(markings[matrix.rowHdr(i).lambdaMarking.id].col)]=true;
	}
	for(std::vector<Marking>::iterator markIter=markingsToKeep.begin() ; markIter!=markingsToKeep.end() ; markIter++)
	{
		assert(checkMarkingValid(*markIter));
		keepCol[matrix.getColNumber(markings[markIter->id].col)]=true;
	}

	// create the list of the columns which should be copied
	std::list<Col> colsToKeep;
	int tmpNumTreedCols=0;
	for(int j=0 ; j<numCols; j++)
	{
		if(keepCol[j])
		{
			colsToKeep.push_back(matrix.getCol(j));
			if(j<(int)numTreedCols)
				tmpNumTreedCols++;
		}
	}

	result->numTreedNodes=numTreedNodes;
	result->numTreedCols=tmpNumTreedCols;
	result->matrix.clone(matrix, colsToKeep);

	//create the mapping, for the new column-handles
	std::vector<Col> oldColNewColMap(numCols);
	for(int j=0 ; j<numCols; j++)
	{
		if(keepCol[j])
		{
			oldColNewColMap[j]=colsToKeep.front();
			colsToKeep.pop_front();
		}
	}

	//create the new nodes
	for(int i=0;i<numRows;i++)
	{
		Marking mark=result->allocateNewMarking(oldColNewColMap[matrix.getColNumber(markings[matrix.rowHdr(i).lambdaMarking.id].col)]);
		Row node=result->matrix.getRow(i);

		int nodeNum;
		// insert node into static list of all nodes
		if(firstDeletedNode == -1)
		{
			nodes.resize(totalNumNodes + 1);
			nodeNum = totalNumNodes;
		}
		else
		{
			nodeNum = firstDeletedNode;
			firstDeletedNode = nodes[nodeNum].nextDeleted;
		}

		result->matrix.rowHdr(node).lambdaMarking=mark;
		result->matrix.rowHdr(node).BV=matrix.rowHdr(i).BV;
		result->matrix.rowHdr(node).indexInNodes=nodeNum;

		markings[mark.id].type=LAMBDA;
		markings[mark.id].node=node;

		nodes[nodeNum].deleted=false;
		nodes[nodeNum].row=node;
		nodes[nodeNum].pc=result;
		totalNumNodes++;
	}

	//allocate the copied markings
	for(std::vector<Marking>::iterator markIter=markingsToKeep.begin() ; markIter != markingsToKeep.end(); markIter++)
	{
		Marking m=result->allocateNewMarking(oldColNewColMap[matrix.getColNumber(markings[markIter->id].col)]);
		*markIter=m;
	}
	return result;
 }


/**
 * Returns the amount of non-zero matrix-enties divided by the total amount of matrix-entries
 */
double PowerCircuitCompMatrix::getMatrixUsage()
{
	int num=0;
	for(unsigned int j=0;j<matrix.getNumCols();j++)
	{
		for(unsigned int i=0;i<matrix.getNumRows();i++)
		{
			if(matrix(i,j)!=ZERO)
				num++;
		}
	}
	return (double)num /((double)matrix.getNumRows()*matrix.getNumCols());
}


/**
 * Prints some statistical data.
 */
void PowerCircuitCompMatrix::printStatistics(std::ostream& os)
{
	os<<"There are currently "<<totalNumMarkings<<" markings and "<<totalNumNodes<<" nodes allocated in all PCs together."<<std::endl;
	os<<"In this PC there are "<<getNumMarkings()<<" markings and "<<getNumNodes()<<" nodes."<<std::endl;
	os<<"Size of the matrix: "<<matrix.getNumRows()<<" x "<<matrix.getNumCols()<<std::endl;
	os<<"Percentage of nonZero entries "<<getMatrixUsage()*100<<". "<<std::endl;
	os<<"Edges in underlying graph: "<<getNumEdges()<<". "<<std::endl<<std::endl;
}

}
