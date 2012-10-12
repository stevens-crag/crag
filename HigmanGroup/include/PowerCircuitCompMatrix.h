/*
 * PowerCircuitCompMatrix.h
 *
 *  Created on: 23.03.2011
 *      Author: ArminW
 */

#ifndef POWERCIRCUITCOMPMATRIX_H_
#define POWERCIRCUITCOMPMATRIX_H_

namespace PC
{

class PowerCircuitCompMatrix : public PowerCircuit
{
private:
	//-------------------------------------------
	// Internal representation of the pc
	struct RowHdr
	{
		Marking 		lambdaMarking;
		bool 			BV;
		unsigned int 	indexInNodes;

		RowHdr() {}
		RowHdr(const Marking& l, bool bv, int i): lambdaMarking(l), BV(bv), indexInNodes(i) {}
	};

	struct ColHdr
	{
		std::list<Marking>		markings;
	};

	// These members hold the actual data of the pc instance
	SignMatrix<RowHdr,ColHdr,unsigned int> matrix;
	unsigned int 			numTreedCols;
	unsigned int 			numTreedNodes;

	//-------------------------------------------
	// Administration of nodes and markings
	// (they are maintained in static arrays)
	enum MarkingType
	{
		NORMAL,LAMBDA
	};

	struct IntNode
	{
		Row						row;
		PowerCircuitCompMatrix* pc;
		bool					deleted;
		int						nextDeleted;
	};

	struct IntMarking
	{
		Col 					col;
		PowerCircuitCompMatrix* pc;
		MarkingType 			type;
		Row 					node; 	// node of which this is the successor marking - only used when type==LAMBDA
		bool 					deleted;
		int						nextDeleted;
	};

	typedef std::vector<IntMarking> MarkingVector;
	typedef std::vector<IntNode> NodeVector;

	static MarkingVector 	markings;
	static NodeVector 		nodes;
	static unsigned int 	totalNumMarkings;
	static unsigned int 	totalNumNodes;
	static int				firstDeletedMarking;
	static int				firstDeletedNode;

	//-------------------------------------------
	// Internal helper methods
	bool 			checkMarkingValid(const Marking& m) const;
	bool 			checkNodeValid(Node n) const;

	void 			deleteNode(Row node);
	void 			deleteCol(Col col);
	void 			deleteMarking(Marking& mark);

	Marking			allocateNewMarking(Col col);
	Row 			newNodeFromMarking(const Marking& mark);

	Row				newOneNode();
	Row 			cloneNode(Row node);

	Marking  		newOneMarking();
	Marking  		newZeroMarking();
	Marking 		newUnitMarking(int onePos);
	Marking 		newCopyColMarking(const Marking& mark);

	void 	 		separateMarkingFromCol(const Marking& m);
	void 			mergeCols(Col targetCol, Col col2);
	int 			getTreedNodePosToGivenCol(Col col);
	void 			moveNodeIntoTreedPartOfMatrix(Row node, int newPos);

	int 			findNewPosOfCompactCol(Col col, bool& equal);

	void 			compactifyFromBottom(Col col);
	void 			insertNewPowerOfTwoNode(unsigned int power);
	unsigned int	calculateCompactRepresentation(int n, Sign result[32]);

	void 			setBV(Row node);
	void 			newDoubleNode(Row node);
	void 			compactify(Col col);

	void 			removeDoubleNodesFromMarkings(Row oldNode, Row newNode);
	int 			insertCompactMarkingIntoTreed(const Marking& mark) ;
	int 			insertCompactColIntoTreed(Col col) ;

	int 			insertNodeIntoTreed(Row node);
	void 			moveColsToTreed(std::list<Col>& colList);
	void 			topSortNode(unsigned int nodeIndex, std::vector<Row>& nodesToSort, std::vector<bool>& visited, std::vector<Row>& topSortPerm);

	void 			extendTree(std::vector<Row>& nodeList, std::vector<Col>& colList);

	void			checkCyclesRecursive(Row node, std::vector<bool> visited);
	void			checkCycles();


	//-------------------------------------------
	// Virtual methods from PowerCircuit, implemented in this class

	virtual void 		incMarkingRefCount(const Marking& mark);
	virtual void 		decMarkingRefCount(Marking& mark);

	virtual	Marking		addMarkings(const Marking& m1, const Marking& m2) ;
	virtual	Marking		invMarking(const Marking& m);
	virtual	Marking		incMarking(const Marking& mark);
	virtual	Marking		intersectMarkings(const Marking& m1, const Marking& m2);
	virtual Marking		cloneMarking(const Marking& mark);
	virtual	Node		cloneNode(Node n);
	virtual Marking		copyMarking(const Marking& m);
	virtual	bool		isMarkingReduced(const Marking& m) const;
	virtual	int			compareMarkings(const Marking& m1, const Marking& m2);
	virtual	int			getRedNodeOrd(Node n);
	virtual int			getNodeSignInMarking(Node n, const Marking& m) const;
	virtual bool		isSuccessorMarking(const Marking& m) const;
	virtual	Node		getIncidentNode(const Marking& m);
	virtual Node		getSmallestNodeInMarking(const Marking& m);
	virtual	Marking		getSuccMarking(Node n);

public:
	// Construction and destruction
							PowerCircuitCompMatrix(int numInitNodes=0, int numInitCols=0);
	virtual		 			~PowerCircuitCompMatrix();

	// Operations on power circuit (override virtual methods of class PowerCircuit)
	virtual	void			reduce();
	virtual	void			reduce(std::vector<Node>& nodeVector, std::vector<Marking>& markingVector);
	virtual	Marking			createMarking(int i = 0);
	virtual	Marking			createMarking(const std::list<Node>& nodes);
	virtual	Marking			createMarking(const std::list<Node>& nodeList, const std::list<int>& signList);
	virtual	Node			createNode(const Marking& succ) ;
	virtual	void			connect(const Marking& m, const Marking& p) ;
	virtual	void 			connectInv(const Marking& m, const Marking& q);
	virtual	void 			remove(const Marking& m) ;
	virtual	Node			getReducedNode(unsigned int ord);
	virtual std::list<Node>	getNodes();
	virtual std::list<Marking>	getMarkings();
	virtual std::list<Node>	getMarkingNodes(const Marking& m);
	virtual	PowerCircuit* 	clone(std::vector<Marking>& markingsToKeep);

	// Statistical data of this power circuit
	virtual void			print(std::ostream& os= std::cout);
	double 					getMatrixUsage();
	virtual void 			printStatistics(std::ostream& os = std::cout);
};

}

#endif /* POWERCIRCUITCOMPMATRIX_H_ */
