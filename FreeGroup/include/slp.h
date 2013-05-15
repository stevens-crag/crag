/**
 * \file slp.h
 * \brief File which includes all definitions for slp module
 *
 * Module implementing algorithms on straight line programs. See \ref slp_description for details.
 */

#pragma once
#ifndef CRAG_FREEGROUP_SLP_H_
#define CRAG_FREEGROUP_SLP_H_

#include "slp_vertex.h"
#include "slp_inspector.h"
#include "slp_vertex_word.h"
#include "arithmetic_sequence.h"
#include "slp_pattern_matching.h"
#include "slp_common_prefix.h"
#include "slp_reduce.h"
#include "slp_vertex_hash.h"
#include "permutation16.h"

/**
 * \page slp_description Straight line programs module description
 *
 * This module implements several algorithms working on so-called Straight Line Programs, or SLPs.
 * \tableofcontents
 *
 * \section vertex Basic structure
 *
 * Each SLP is represented as a set of vertices with some root. Each vertex is an object of type Vertex, and
 * there are three implementations of this interface - empty vertex (producing empty word), terminal vertex
 * (producing just one letter) and non-terminal vertex (having left and right child and concatenating
 * words produced by children. Interface documentation: \link crag::slp::Vertex\endlink
 *
 * \subsection examples Usage examples
 * \subsubsection null_vertex Null vertex
 * Vertex() represents empty word. Basically you never need to use it, here are its properties:
 * \code
 *   Vertex().length() == 0;
 *   Vertex().height() == 0;
 *   Vertex().left_child() == Vertex();
 *   Vertex().negate() == Vertex();
 * \endcode
 *
 * \subsubsection terminal_vertex Terminal vertex
 * Terminal Vertices are implemented as a class with template parameter TerminalSymbol, so one can vary
 * the group alphabet. Usually only one alphabet is used, so we recommend to typedef it in your
 * functions:
 * \code
 *   typedef TerminalVertexTemplate<char> TerminalVertex; //We are going to use char as group alphabet
 *   typedef TerminalVertexTemplate<int> TerminalVertexOnInfiniteAlphabet; //We also could use int as a base type, or anything else, which meets certain criteria
 * \endcode
 *
 * To get a vertex for a terminal symbol, just create the object of its type:
 * \code
 *   TerminalVertex a('a'); //a represents $a$
 * \endcode
 *
 * Properties of TerminalVertexTemplate<TerminalSymbol>:
 * \code
 *   a.length() == 1; //TRUE
 *   a.height() == 1; //TRUE
 *   a.left_child() == Vertex(); // a has no children, so return empty vertex
 * \endcode
 *
 * From object of type TerminalVertex you also can get the underlying alphabet symbol:
 * \code
 *   a.terminal_symbol() == 'a'; //you can use member function which will return you stored symbol
 *   static_cast<char>(a) == 'a'; //also TerminalVertex can be casted to its basic type, but only explicitly
 * \endcode
 *
 * Your often have to get the vertex representing the inverse, i.e. \f$w^{-1}\f$. Terminal Vertex just return another TerminalVertex which store
 * -terminal_symbol():
 * \code
 *   TerminalVertex a1 = a.negate(); //a1 represents $a^{-1}$.
 *   a1.terminal_symbol() == -'a';
 * \endcode
 *
 * The equality operator is defined for all vertices. However, vertices of different types (i.e., Vertex::Null and TerminalVertex) are not equal,
 * even if they producing the same words. Inside each type they are compared in a bit different way. TerminalVertex are equal if and only if they store
 * same terminal symbols:
 * \code
 *   a == a;
 *   a.negate().negate() == a;
 *   a.negate() != a;
 *   TerminalVertex('a') == TerminalVertex('a');
 *   TerminalVertex('a') != TerminalVertex('b');
 * \endcode
 *
 * \subsubsection nonterminal_vertex Non-terminal vertex
 *
 * To construct some non-terminal vertex, create an object of NonterminalVertex:
 * \code
 *   TerminalVertex b('a'); //b produces $b$
 *   NontermianlVertex ab(a, b); //This vertex represents word $ab$
 *   NonterminalVertex aa(a, a); //You can use the same vertices as a children
 *   NonterminalVertex aba1a1(ab, aa.negate()); //Here is something complex
 *   NonterminalVertex(a, a.negate()) != Vertex::Null; //Words are not reduces automatically
 * \endcode
 *
 * Properties of NonterminalVertex (this time they depends on the children):
 * \code
 *   ab.length() == 2;
 *   ab.height() == 2;
 *   ab.left_child() == a;
 *   ab.right_child() == b;
 *
 *   aba1a1.length() == 4;
 *   aba1a1.height() == 3;
 *   aba1a1.left_child() == ab;
 *   aba1a1.right_child() == a1a1.negate();
 *   aba1a1.right_child().left_child() == a1;
 * \endcode
 *
 * As you already have seen, non-terminal vertices can also be 'negated':
 * \code
 *   NontermianlVertex b1a1 = ab.negate(); //represents $(ab)^{-1} = b^{-1}a^{-1}$
 *   b1a1.left_child() == b.negate();
 *   b1a1.right_child() == a1;
 *   b1a1.negate().left_child() == a; // b1a1.negate() is $ab$
 * \endcode
 *
 * When you compare NonteriminalVertices, you do not compare the words, produced by them, but rather the vertices objects. Each time you create a new
 * NonterminalVertex(Vertex, Vertex), a new vertex is created. This is not true for NonterminalVertex(NonterminalVertex) and for negate():
 * \code
 *   ab != NonterminalVertex(a, b);
 *   ab.negate().negate() == ab;
 *   ab.negate() == ab.negate();
 *   ab == NonterminalVertex(ab);
 * \endcode
 *
 * left_child() and right_child() also produce new vertex each time:
 * \code
 *   aba1a1.left_child() != aba1a1.left_child()
 * \endcode
 *
 * \subsubsection general General Vertex
 * All these types just implement the common interface defined in class Vertex. So, you can write some functions working with all these types:
 * \code
 *   //Find the leftmost terminal which is a child of some vertex
 *   Vertex leftmost_terminal_vertex(Vertex vertex) {
 *     while (vertex.height() > 1) { //only terminal vertices has height 1
 *       //vertex == false only if vertex == Vertex()
 *       //if vertex has no left child, go to right
 *       vertex = vertex.left_child() ? vertex.left_child() : vertex.right_child();
 *     }
 *     return vertex;
 *   }
 * \endcode
 *
 * If you want then to get a terminal vertex from this function (for example, to get terminal symbol), you can construct
 * TerminalVertexTemplate<TerminalSymbol>(Vertex):
 * \code
 *   TerminalVertex<char> t = leftmost_terminal_vertex(aba1a1); //t should represent $a$
 *   t.terminal_symbol() == 'a';
 * \endcode
 *
 * If you will try to cast non-terminal vertex to terminal vertex, you will get Vertex::Null:
 * \code
 *   TerminalVertex<char>(aa) == Vertex();
 *   TerminalVertex<char>(Vertex()) == Vertex();
 * \endcode
 *
 * \section vertex_word Words produced by SLP
 *
 * VertexWord<TerminalSymbol> interprets vertex of SLP as a string of TerminalSymbol, therefore allowing you to take specific symbols, or iterate
 * over a word represented by a vertex.
 *
 * \subsection usage Usage
 * You just create VertexWord<TerminalSymbol> giving Vertex as an argument:
 * \code
 *   VertexWord<char> word(aba1a1);
 * \endcode
 * To retrieve specific symbol, use operator[]:
 * \code
 *   word[0] == 'a';
 *   word[3] == -'a';
 * \endcode
 * You also can iterate through the word:
 * \code
 *   auto i = word.begin();
 *   *i == 'a';
 *   ++i;
 *   *i == 'b';
 * \endcode
 * The iterator can be moved only forward, but it supports move for more than one symbol at time:
 * \code
 *   i += 2'
 *   *i == word[3];
 *   auto j = word.begin() + 2;
 *   *j == word[2];
 * \endcode
 *
 * It also can be used as standard iterators in most cases:
 * \code
 *   for (auto k = word.begin(); k != word.end(); ++k) {
 *     ::std::cout << *k;
 *   }
 *
 *   for (auto symbol : word) { //new 'foreach' semantics
 *     ::std::cout << symbol;
 *   }
 *
 *   ::std::string word_string(word.begin(), word.end());
 *   word_string[0] == 'a';
 *   word_string[1] == 'b';
 *   word_string[2] == -97; //-'a'
 *   word_string[3] == -'a';
 * \endcode
 *
 * The equality operator is defined for VertexWord. It do not compare words symbol-by-symbol (which can take exponential time), but
 * use pattern_matching.h algorithms.
 *
 * \section inspector SLP Traversal
 *
 * All tree traversal routing is implemented in slp_inspector.h module. The basic usage is following:
 * \code
 *   using crag::slp;
 *
 *   PostorderInspector inspector(slp_vertex);
 *
 *   while (!inspector.stopped()) {
 *     //some code which can refer to the current vertex using inspector.vertex()
 *     ...
 *     inspector.next()
 *   }
 * \endcode
 *
 * You have several possibilities to customize its behavior. First of all, you can define traversal order by
 * using one of crag::slp::PreorderInspector, crag::slp::InorderInspector, crag::slp::PostorderInspector. Second,
 * you can 'omit' some vertices considering them as 'Null' vertices (e.g. all vertices of height < 10). Note that you
 * also will skip all their children. As an example, we provide the following very useful code. It visit every vertex
 * only once, even if it is a child of more then one vertex:
 * \code
 *   using crag::slp;
 *
 *   std::unordered_set<Vertex> visited_vertices;
 *   auto acceptor = [&visited_vertices] (const inspector::InspectorTask& task) {
 *     return visited_vertices.count(task.vertex) == 0; //true only if vertex is not visited yet
 *   };
 *   Inspector<inspector::Postorder, decltype(acceptor)> inspector(root, acceptor);
 *   while (!inspector.stopped()) {
 *     visited_vertices.insert(inspector.vertex());
 *     //some code which deals with inspector.vertex()
 *     ...
 *     inspector.next()
 *   }
 * \endcode
 *
 * Further customization can be done using #crag::slp::inspector::InspectorPath policy. See source code for examples.
 *
 * \section matchin_table Pattern matching using matching tables
 *
 * Submodule slp_pattern_matching.h implements the pattern matching algorithms based on crag::slp::MatchingTable.
 * Matching table itself have the following interface:
 * \code
 *   using crag::slp;
 *   MatchingTable table;
 *   crag::FiniteArithmeticSequence match = table.matches(pattern, text); //pattern, text are of type Vertex
 * \endcode
 *
 * \a match now keeps the arithmetic sequence of starting position of pattern inside text, which touch text.split_point().
 *
 * If you want to inspect all inclusions of pattern into text, you will be served with crag::slp::PatternMatchesGenerator.
 * Here is usage:
 *
 * \code
 *   crag::slp::PatternMatchesGenerator matches(pattern, text);
 *   do {
 *     auto next_match = matches.next_match();
 *     //do something with next_match
 *   } while (matches);
 *   //post check is more convenient here, since matches == true also just after last match was generated
 * \endcode
 *
 * \section misc Small(and a bit larger) miscellaneous features
 *
 * \subsection Getting sub slp
 *
 * One useful function is defined in slp_reduce.h which is call crag::slp::get_sub_slp(const Vertex& root, const LongInteger& begin, const LongInteger& end).
 * This function carefully (i.e. thinking befor copying) construct an slp which generate a subword starting at begin and ending before end.
 *
 * \subsection Longest common prefix
 *
 * If you need to find the longest common prefix of two sequences, you can use crag::slp::longest_common_prefix(const Vertex& first, const Vertex& second)
 * from slp_common_prefix.h, which returns the though-for value.
 *
 * \subsection Reduction
 *
 * If you want to cancel out all substrings like \f$cc^{-1}\f$, call crag::slp::reduce(const Vertex& root) from slp_reduce.h
 *
 * \subsection Fast longest common prefix & reduction with hashing
 *
 * If speed of metioned methods will be too slow for you, we also have implemented some algorithms in slp_vertex_hash.h. They are
 * not documented yet, ask developers for its usage.
 */
#endif /* CRAG_FREEGROUP_SLP_H_ */
