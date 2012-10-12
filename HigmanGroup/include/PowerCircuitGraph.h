/*
 * PowerCircuitCompMatrix.h
 *
 *  Created on: 23.03.2011
 *      Author: ArminW
 */

#ifndef POWERCIRCUITGRAPH_H_
#define POWERCIRCUITGRAPH_H_


namespace PC
{

class PowerCircuitGraph : public PowerCircuit
{
private:


	//-------------------------------------------
	// Administration of nodes and markings

	struct IntNode
	{
		bool 					BV;
		unsigned int 			indexInOrder;
		bool					deleted;
		int						nextDeleted;
		std::list< std::pair<Node,Sign> > 	successors;
	};

	struct IntMarking
	{
		bool 					deleted;
		bool					reduced;
		int						nextDeleted;
		int						refCount;
		std::list< std::pair<Node,Sign> > 	nodes;
	};

	struct NodeUsedByType
	{
		std::list< std::pair<Node,Sign> >* 			nodeList;
		std::list< std::pair<Node,Sign> >::iterator	position;
		int 										id;
		NodeUsedByType(std::list< std::pair<Node,Sign> >* l, std::list< std::pair<Node,Sign> >::iterator iter, int i )
		:nodeList(l), position(iter), id(i){}

		bool operator== (const NodeUsedByType& t2)
		{
			return (id==t2.id);
		}
	};


	// These members hold the actual data of the pc instance
	unsigned int 			numNodes;
	unsigned int 			numMarkings;
	unsigned int 			numReducedNodes;
	std::vector<Node>		nodeOrder;

	std::vector<IntNode>	nodes;
	std::vector<IntMarking>	markings;

	int						firstDeletedNode;
	int						firstDeletedMarking;


	/*
	 * static method, to make the sort algorithm on signed lists run
	 */
	static bool	compareNodesLessThan(const std::pair<Node,Sign>& n1, const std::pair<Node,Sign>& n2)
	{
		assert(n1.first.pc==n2.first.pc);
		PowerCircuitGraph* pc = (PowerCircuitGraph*) n1.first.pc;
		return (pc->nodes[n1.first.id].indexInOrder < pc->nodes[n2.first.id].indexInOrder);
	}


	//-------------------------------------------
	// Internal helper methods
	bool 			checkMarkingValid(const Marking& m) const;
	bool 			checkNodeValid(Node n) const;

	void 			deleteNode(Node node);
	void 			deleteMarking(Marking&  mark);

	Marking			newMarking(std::list< std::pair<Node,Sign> >  nodeList );
	Node 			newNode(std::list< std::pair<Node,Sign> >  nodeList );

	Node			newOneNode();

	Marking  		newOneMarking();
	Marking  		newZeroMarking();
	Marking 		newUnitMarking(int onePos);
	Marking 		newCopyMarking(const Marking& mark);

	std::list< std::pair<Node,Sign> >	inc(const std::list< std::pair<Node,Sign> >& nodeList);

	void 			moveNodeIntoReducedPart(Node node, int newPos);

	int 			compare(std::list< std::pair<Node,Sign> >& l1, std::list< std::pair<Node,Sign> >&  l2);
	void 			sortNodeList(std::list< std::pair<Node,Sign> >&  l);

	void 			insertNewPowerOfTwoNode(unsigned int power);
	std::list< std::pair<Node,Sign> >	calculateCompactRepresentation(int n);

	void 			markSucessors(std::vector<bool> marked, std::list< std::pair<Node,Sign> >& nodeList);
	void 			setBV(Node node);
	Node 			newDoubleNode(Node node);

	std::list< std::pair<Node,Sign> >::iterator	findNodeInList(unsigned int nodeIndex, std::list< std::pair<Node,Sign> >& nodeList);
	void 			removeDoubleNodesFromMarkings(Node oldNode, Node newNode, std::list<NodeUsedByType>* nodeUsedBy);

	int 			findNewPosOfReducedNode(Node node, bool& equal);
	int 			insertNodeIntoReduced(Node node, std::list<NodeUsedByType>* nodeUsedBy);
	void 			topSortNode(Node node, std::vector<bool>& visited, std::vector<Node>& topSortPerm);

	void 			extendTree(std::vector<Node>& nodeList, std::vector<Marking>& markingList);

	//-------------------------------------------
	// Virtual methods from PowerCircuit, implemented in this class
	virtual void 		incMarkingRefCount(const Marking& mark);
	virtual void 		decMarkingRefCount(Marking& mark);
	virtual	Marking		addMarkings(const Marking& m1, const Marking& m2) ;
	virtual	Marking		invMarking(const Marking&  m);
	virtual	Marking		incMarking(const Marking&  mark);
	virtual	Marking		intersectMarkings(const Marking&  m1, const Marking&  m2);
	virtual Marking		cloneMarking(const Marking&  mark);
	virtual	Node		cloneNode(Node n);
	virtual Marking		copyMarking(const Marking&  m);
	virtual bool 		isMarkingReduced(const Marking&  m) const;
	virtual	int			compareMarkings(const Marking&  m1, const Marking&  m2);
	virtual	int			getRedNodeOrd(Node n);
	virtual int			getNodeSignInMarking(Node n, const Marking&  m) const;
	virtual bool		isSuccessorMarking(const Marking&  m) const;
	virtual	Node		getIncidentNode(const Marking&  m);
	virtual Node		getSmallestNodeInMarking(const Marking&  m);
	virtual	Marking		getSuccMarking(Node n);

public:
	// Construction and destruction
							PowerCircuitGraph();
	virtual		 			~PowerCircuitGraph();

	// Operations on power circuit (override virtual methods of class PowerCircuit)
	virtual	void			reduce();
	void 					reduce(std::list<Marking> markingList);
	virtual	void			reduce(std::vector<Node>& nodeVector, std::vector<Marking>& markingVector);
	virtual	Marking			createMarking(int i = 0);
	virtual	Marking			createMarking(const std::list<Node>& nodes);
	virtual	Marking			createMarking(const std::list<Node>& nodeList, const std::list<int>& signList);
	virtual	Marking			createMarking(const std::list< std::pair<Node,Sign> >& nodeList );
	virtual	Node			createNode(const Marking& succ) ;
	virtual	void			connect(const Marking& m, const Marking&  p) ;
	virtual	void 			connectInv(const Marking& m, const Marking&  q);
	virtual	void 			remove(const Marking& m) ;
	virtual	Node			getReducedNode(unsigned int ord);
	virtual std::list<Node>	getNodes();
	virtual std::list<Marking>	getMarkings();
	virtual std::list<Node>	getMarkingNodes(const Marking& m);
	virtual	PowerCircuit* 	clone(std::vector<Marking>& markingsToKeep);

	// Statistical data of this power circuit
	virtual void			print(std::ostream& os= std::cout);
	void  					printMarking(const Marking& m, std::ostream& os= std::cout);
	void 					printMarking(const std::list< std::pair<Node, Sign> >& l ,std::ostream& os= std::cout);
	double 					getMatrixUsage();
	virtual void 			printStatistics(std::ostream& os = std::cout);

	bool					checkCyclesRecursive(int n, std::vector<bool> & visited);
	bool					checkCycles();
	bool 					checkConsistency();
};

}

#endif /* POWERCIRCUITGRAPH_H_ */
