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
 * Vertex::Null represents empty word. Basically you never need to use it, here are its properties:
 * \code
 *   TerminalVertex::Null.length() == 0;
 *   TerminalVertex::Null.height() == 0;
 *   TerminalVertex::Null.left_child() == TerminalVertex::Null;
 *   TerminalVertex::Null.negate() == TerminalVertex::Null;
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
 *   a.left_child() == Vertex::Null; // a has no children, so return empty vertex
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
 *       //vertex == false only if vertex == Vertex::Null
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
 *   TerminalVertex<char>(aa) == Vertex::Null;
 *   TerminalVertex<char>(Vertex::Null) == Vertex::Null;
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
 */

#endif /* CRAG_FREEGROUP_SLP_H_ */
