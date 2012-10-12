/*
 * PowerCircuit.h
 *
 *  Created on: 23.03.2011
 *      Authors: Juern, ArminW
 */


/**
 * 
 *
 * This package implements power circuits in two different ways. The first one
 * uses compact representation of markings as defined in [1]. The tree containing
 * all markings is stored not as a graph but as a matrix to improve access time. 
 * This implementation is asymptotically optimal. 
 *
 * The other implementation (PowerCircuitGraph) stores the graph by adjacency lists
 * and resembles the less technical approach chosen in [2].
 * Markings, too, are stored as lists. This implementation does not make use of compact
 * representations. Thus, it is not theoretically optimal, but usually faster for
 * typical inputs as it uses a lot of overhead.
 *
 * Both power-circuit-classes are derived from the abstract class PowerCircuit and
 * can be used interchangably.
 *
 * To perform calculations with power circuits, create an instance of one of the above
 * classes. The PowerCircuit class provides methods to create markings, perform
 * reduction, insert edges into the graph (connect), and other tasks.
 *
 * Markings work like references to the respective internal representations. Thus the
 * assignment operator does not create a deep copy of a marking. However all arithmetic
 * operations create new markings, so you can keep the old ones. It is also not necessary
 * (and not possible) to delete markings (in the PowerCircuitGraph class once there is
 * no marking pointing onto an internal marking, its memory is automatically deallocated).
 * Just make sure that you create markings only for the time you really need them.
 *
 * To get an idea of how to use power circuits, you might want to take a look at the
 * examples included in this package, where power circuits are used to solve the word
 * problem in certain groups. 
 *
 * In order to use the draw method you need graphviz installed on your computer (download from
 * http://www.graphviz.org/). This method creates a postscript file with a graphical depiction
 * of the graph.
 *
 * The documentation is created with doxygen (see http://www.stack.nl/~dimitri/doxygen/index.html).
 * To create the documentation in html and latex format use "doxygen".
 * You can change the settings for the documentation in the file "Doxyfile". If "Doxyfile"
 * does not exist, use "doxygen -g" to create it.
 *
 * \authors Juern Laun
 * \authors Armin Weiss
 *
 * A lot of the theoretical background of power circuits is due to Sasha Ushakov. We thank him
 * for his helpful advice, and also for the invitation to th Stevens Institute where parts of this
 * implementation were developed. 
 *
 * [1] V. Diekert, J. Laun, A. Ushakov. Efficient algorithms for highly compressed data: The Word
 * Problem in Higman's Group is in P. Preprint, 2011. http://arxiv.org/abs/1103.1232
 *
 * [2] V. Diekert, J. Laun, A. Ushakov. Efficient algorithms for highly compressed data: The Word
 * Problem in Higman's Group is in P. Conference version, to appear, 2011.
 */

#ifndef POWERCIRCUIT_H_
#define POWERCIRCUIT_H_

namespace PC
{

/*
 * Forward declarations of the main classes
 */
class PowerCircuit;
class Marking;
class Node;

/**
 * Wrapper class for nodes in power circuits. In
 * reality, nodes are handled by (derivatives of)
 * the class PowerCircuit. This class is used to
 * make handling of nodes type-safe. It contains
 * only a pointer to the corresponding power
 * circuit and an integer id, which is -1 if the
 * node has not (yet) been set. If any more data
 * is needed for a node, the power circuit class
 * has to hold it itself. The methods call
 * similarly-named methods of PowerCircuit.
 */
class Node
{
public:
	PowerCircuit*		pc;
	int					id;

						Node(): pc(NULL), id(-1) {}
						Node(PowerCircuit* p, int i): pc(p), id(i) {}
	bool				isReduced();
	bool				isDefined();
	Node				clone();
	int					getOrder();
	Marking				getSuccessorMarking();
	bool 				operator==(Node n);

};

/**
 * Wrapper class for markings in power circuits. In
 * reality, markings are handled by (derivatives of)
 * the class PowerCircuit. This class is used to
 * make handling of markings type-safe. It contains
 * only a pointer to the corresponding power
 * circuit and an integer id, which is -1 if the
 * marking has not (yet) been set. If any more data
 * is needed for a marking, the power circuit class
 * has to hold it itself. The methods call
 * similarly-named methods of PowerCircuit.
 *
 * Markings work like references to the respective internal
 * representations. Thus the assignment operator does not
 * create a copy of a marking. However all arithmetic
 * operations create new markings, so the old ones can be maintained.
 * Therefore a marking (i.e. the handle of a internal representation)
 * only changes its value through the assignment operator and
 * when calling the connect() or connectInv() methods of
 * PowerCircuit.
 *
 * It is also not necessary (and not possible) to delete markings
 * (in the PowerCircuitGraph class, once there is no marking pointing
 * onto an internal marking, its memory is automatically deallocated).
 *
 */
class Marking
{
public:
	PowerCircuit*		pc;
	int					id;

						Marking(): pc(NULL), id(-1) {}
						Marking(PowerCircuit* p, int i);
						Marking(const Marking& m);
						Marking& operator= (const Marking& rhs);
						~Marking();

	bool				isReduced();
	Marking				clone();
	Marking 			copy();
	std::list<Node>		getNodes();
	int					getSign(Node n);
	Marking				operator++(int dummy);
	Marking				operator-();
	Marking				operator+(const Marking& m);
	Marking				operator&(const Marking& m);
	bool 				operator<(const Marking& m);
	bool 				operator>(const Marking& m);
	bool 				operator<=(const Marking& m);
	bool 				operator>=(const Marking& m);
	bool 				operator==(const Marking& m);
	bool 				operator!=(const Marking& m);
	bool				isSuccessorMarking();
	bool				isDefined();
	Node				getIncidentNode();
	Node				getSmallestNode();
};

/**
 * Abstract base class for all implementations of power circuits.
 * The PowerCircuit class provides methods to create markings, perform the reduction procedure,
 * draw edges in the graph (connect), and some additional tasks.
 */
class PowerCircuit
{
private:
	// Internal helper methods (these are mostly called by the classes Node and Marking)

	virtual void 			incMarkingRefCount(const Marking& mark) = 0;
	virtual void 			decMarkingRefCount(Marking& mark) = 0;

	virtual	Marking			addMarkings(const Marking& m1, const Marking&  m2) = 0;
	virtual	Marking			invMarking(const Marking& m) = 0;
	virtual	Marking			incMarking(const Marking& m) = 0;
	virtual	Marking			intersectMarkings(const Marking& m1, const Marking&  m2) = 0;
	virtual Marking			cloneMarking(const Marking& m) = 0;
	virtual	Node			cloneNode(Node n) = 0;
	virtual Marking			copyMarking(const Marking& m) = 0;
	virtual	bool			isMarkingReduced(const Marking& m) const = 0;
	virtual	int				compareMarkings(const Marking& m1, const Marking&  m2) = 0;
	virtual	int				getRedNodeOrd(Node n) = 0;
	virtual int				getNodeSignInMarking(Node n, const Marking&  m) const = 0;
	virtual bool			isSuccessorMarking(const Marking& m) const = 0;
	virtual	Node			getIncidentNode(const Marking& m) = 0;
	virtual Node			getSmallestNodeInMarking(const Marking& m) = 0;
	virtual	Marking			getSuccMarking(Node n) = 0;


	friend class Node;
	friend class Marking;

public:
	// Construction and destruction
							PowerCircuit() {}
	virtual					~PowerCircuit() {}


	// Operations on power circuit

	/**
	 * Reduces the whole PC and all its markings.
	 */
	virtual	void			reduce() = 0;

	/**
	 *  Moves all nodes in nodeList and markings in markingList to the reduced part of the PC.
	 *  Nodes in nodeList MUST NOT be reduced!
	 *  The support of markings in markingList has to be in nodeList or the already reduced part of the PC.
	 */
	virtual	void			reduce(std::vector<Node>& nodeVector, std::vector<Marking>& markingVector) = 0;

	/**
	 * Creates a new marking with value val. The therefore required nodes are inserted into the reduced part of the PC.
	 */
	virtual	Marking			createMarking(int i = 0) = 0;

	/**
	 *  Creates a new marking, with the nodes in nodeList set to 1.
	 */
	virtual	Marking			createMarking(const std::list<Node>& nodes) = 0;

	/**
	 *  Creates a new marking, with the nodes in nodeList set to the respective sign in signList.
	 */
	virtual	Marking			createMarking(const std::list<Node>& nodeList, const std::list<int>& signList) = 0;

	/**
	 *  Creates a new marking, with the nodes set to 1. As first argument it requires the number of nodes set to one.
	 *  The other arguments are the respective nodes.
	 */
	virtual Marking			createMarkingFromNodes(unsigned int numNodes, ...);

	/**
	 * Creates a new node with a copy of m as successor-marking.
	 */
	virtual	Node			createNode(const Marking& succ) = 0;

	/**
	 * Adds the marking p to all successor markings of nodes in the support of m. (So, if supp(m) and supp(p) are
	 * disjoint, that means, that arrows are drawn from all nodes in supp(m) to p).
	 * The user has to assure that, for each node, the sum lies in the range [-1,+1] (for example using clone())!
	 * The user is required to check that the resulting graph has no cycles and no successor marking becomes negative.
	 * All nodes in supp(m) have to be non-reduced (use clone() to assure that).
	 */
	virtual	void			connect(const Marking& m, const Marking&  p) = 0;

	/**
	 * The same as connect, only that the marking p is subtracted instead of added.
	 */
	virtual	void			connectInv(const Marking& m, const Marking&  p) = 0;

	/**
	 * Deletes all nodes in the given marking. (So if other markings use those
	 * nodes, their values may be changed!)
	 */
	virtual void			remove(const Marking& m) = 0;

	/**
	 * Returns the node with index ord (that means there are ord smaller nodes in the reduced part of the PC)
	 * If there are less than ord nodes, the undefined node is returned.
	 */
	virtual	Node			getReducedNode(unsigned int ord) = 0;

	/**
	 * Returns a list of all nodes in the PowerCircuit.
	 */
	virtual std::list<Node>	getNodes() = 0;

	/**
	 * Returns a list of all markings in the PowerCircuit.
	 */
	virtual std::list<Marking>	getMarkings() = 0;

	/**
	 * Returns a list of all nodes which are non-zero in the marking m.
	 */
	virtual std::list<Node>	getMarkingNodes(const Marking& m)  = 0;

	 /**
	  * Makes a copy of a PowerCircuit. Only the markings in the markingsToClone list are available in the cloned circuit.
	  * Reduced nodes and markings are remain reduced.
	  * The clone has to be deallocated via delete, once it is not needed anymore!
	  */
	virtual PowerCircuit*	clone(std::vector<Marking>& markingsToClone) = 0;

	virtual	void			draw(std::string filename, Marking highlight1 = Marking(), Marking highlight2 = Marking(), Marking highlight3 = Marking(), Marking highlight4 = Marking(), Marking highlight5 = Marking(), Marking highlight6 = Marking(), Marking highlight7 = Marking(), Marking highlight8 = Marking(), Marking highlight9 = Marking());

	// Statistical data of this power circuit
	int 					getNumEdges();
	int 					getNumNodes();
	int 					getNumMarkings();

	/**
	 * Prints some statistical data.
	 */
	virtual void 			printStatistics(std::ostream& os = std::cout);

	/**
	 * Prints the adjacency matrix of the underlying graph.
	 */
	virtual void			print(std::ostream& os= std::cout) = 0;
};

}

#endif /* POWERCIRCUIT_H_ */

