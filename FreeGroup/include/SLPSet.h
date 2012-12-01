/*
 * File:   CompositionSystemSet.h
 * Author: dpantele
 *
 * Created on November 18, 2012, 4:06 PM
 */

#ifndef COMPOSITIONSYSTEMSET_H
#define	COMPOSITIONSYSTEMSET_H

#include <gmpxx.h>
typedef mpz_class LongInteger;

#include <vector>
#include <memory>
#include <initializer_list>

//! Represents a vertex in program or the "inversed" vertex.
struct SignedVertex {
  size_t index;     //!< The index of the vertex in in the vertices vector.
  bool is_negative; //!< True if we need "inversed" vertex

  //Constructor to simplify the creation of this vertex
  SignedVertex(size_t index, bool negative)
          : index(index)
          , is_negative(negative)
  { }

};

bool operator==(const SignedVertex& lhs, const SignedVertex& rhs) {
  return lhs.index == rhs.index && lhs.is_negative == rhs.is_negative;
}



//! Struct representing one vertex in the SLP. Internal.
/**
 * Represents one composition rule of kind \p$A\rightarrow BC \p$.
 * Also stores internal information such as the height of the subtree
 * and the length of it.
 *
 * TODO: check if we can guarantee that any non-terminal vertex has exactly
 * two children.
 *
 */
struct SLPVertex {
  //! Constant indicating that there is no such child.
  const SignedVertex CHILD_NOT_EXIST = SignedVertex(-1, false); 
  SignedVertex left_child;               //!< The index of the left vertex in vertices vector. Use CHILD_NOT_EXIST
  SignedVertex right_child;              //!< The index of the right vertex in vertices vector.
  unsigned int terminal_symbol;          //!< NON_TERMINAL, if non-terminal. Otherwise the number of the symbol, greater that zero.
  const unsigned int NON_TERMINAL = 0;   //!< Constant, meaning this vertex is non-terminal
  LongInteger length;                    //!< Length of word produced by the vertex
  unsigned int parents_count;            //!< Number of parents

  SLPVertex()
    : left_child(CHILD_NOT_EXIST)
    , right_child(CHILD_NOT_EXIST)
    , terminal_symbol(NON_TERMINAL)
    , length(0)
    , parents_count(0)
  { }
};

//! Progression tables which are described the the thesis by Lifshits
  /**
   * It is a class of the progression table as described in the thesis
   * by Yury Lifshits. This table enumerates the entries of the subtrees
   * of pattern to the subtrees of text.
   *
   * It is calculated for two SLPs \em P and \em T, and has \em nm cells, where
   * \em n is the number of vertices in \em P and \em m is the number of vertices
   * in \em T. In the cell PT[i][j] we have all entries of the word produced by
   * the vertex \p$P_i\p$ into the word produced by vertex \p$T_j\p$, which have
   * some common part with the <em>split point</em> of \p$T_j\p$. Split point is
   * the position in the word \p$T_j\p$ just after the end of the first part,
   * i.e. part \p$T_r\p$ of the production rule \p$T_j \to T_r T_s\p$.
   *
   * The important fact that we use here is that if some entries have the common
   * point (split point in this case), then the beginning of these entries make
   * the arithmetic progression, which can be encoded by the triplet of integers.
   * So, the table just stores these triplets.
   *
   * See the "Algorithms and complexity analysis for processing compressed text"
   * for details.
   *
   * @param pattern 
   * @param text The SLP of the text where we are looking the word.
   * @param prefix_table The pre-allocated storage for the result.
   */
class ProgressionTable {
public:
  struct MatchResultSequence {
    LongInteger start; //!< The beginning of the first match
    LongInteger step;  //!< The distance between the matches
    LongInteger count; //!< The number of matches
  };

  //! Return all matches of the pattern around the "split point"
  /**
   * This function get the result from #table and recursively calculate
   * it if needed.
   * 
   * @param pattern The SLP for a word we want to find in text.
   *                Specified by the index of vertex in SLPSet::vertices.
   *
   * @param text    The SLP for the text where we are searching.
   *                Index in SLPSet::vertices
   * @return The sequence of the beginnings of matches.
   */
  MatchResultSequence get_matches(const SignedVertex& pattern,
                                  const SignedVertex& text);

  ProgressionTable(const std::vector<SLPVertex>& pattern_vertices,
                   const std::vector<SLPVertex>& text_vertices)
      : pattern_vertices(pattern_vertices)
      , text_vertices(text_vertices)
      , table(pattern_vertices.size() * text_vertices.size())
  {}



private:
  struct MatchResult {
    MatchResultSequence match; //!< Entries of pattern in text
    MatchResultSequence inversed_match; //!< Entires of inversed pattern in text
  };

  //! Helper function which looks for pattern in text[begin..end]
  std::pair<MatchResultSequence, MatchResultSequence>
    local_search(SignedVertex pattern, SignedVertex text,
                 LongInteger begin, LongInteger end);


  const std::vector<SLPVertex>& pattern_vertices; //!< Reference to SLP with pattern
  const std::vector<SLPVertex>& text_vertices; //!< Reference to SLP with text
  std::vector<MatchResult> table; //!<The storage for matchings, virtual 2D
};

//! class of straight line program collection
/**
 * A straight line program is a program which computes a single word.
 * We used "Polynomial-time Word problems" by Saul Schleimer and 
 * "Algorithms and complexity analysis for processing compressed text"
 * by Yury Lifshits when implementing this class
 *
 * Actually, this class computes several words at the same time, one word for
 * each "root" of this system.
 */
class SLPSet {
public:
  SLPSet(size_t terminals_count)
    : vertices(terminals_count)
    , roots(terminals_count)
  {
    for (size_t terminal_id = 0; terminal_id < terminals_count; ++terminal_id) {
      vertices[terminal_id].terminal_symbol = terminal_id + 1;
      vertices[terminal_id].length = 1;
      vertices[terminal_id].parents_count = 1;
      roots[terminal_id] = terminal_id;
    }
  }

  //! Constructor to create almost trivial SLP with one non-trivial root
  /**
   * Constructs SLP with #terminals_count terminals and the same number of roots.
   * Root #nontrivial_root (index starts from 1) has the only non-trivial
   * production rule (rule_lhs, rule_rhs).
   * 
   * @param terminals_count The number of terminals and roots
   * @param nontrivial_root The index of non-trivial root, starting from 1
   * @param rule_lhs The terminal index of the left part of the only non-trivial
   *                 production rule. If < 0, then terminal is reversed.
   * @param rule_rhs The right side of prudction rule, see docs for rule_lhs
   */

  SLPSet(size_t terminals_count, size_t nontrivial_root, int rule_lhs, int rule_rhs)
    : SLPSet(terminals_count)
  {
    SLPVertex non_trivial();

    if (rule_lhs != -rule_rhs) {//If one vertex cancels another one, leave just root without any children
      non_trivial.left_child.index  = abs(rule_lhs) - 1;
      non_trivial.right_child.index = abs(rule_rhs) - 1;
      non_trivial.left_child.is_negative = (rule_lhs < 0);
      non_trivial.right_child.is_negative = (rule_rhs < 0);
      non_trivial.length = 2;
      if (abs(rule_lhs) != nontrivial_root) {
        vertices[non_trivial.left_child.index].parents_count += 1;
      }
      if (abs(rule_rhs) != nontrivial_root) {
        vertices[non_trivial.right_child.index].parents_count += 1;
      }
    }
    
    vertices.push_back(non_trivial);
    roots[nontrivial_root - 1] = vertices.size() - 1;

  }

  //! Compose current program with another one
  /**
   * Add all vertices from the other SLP to *this and connect \p$i\p$-th
   * terminal of the other SLP to \p$i\p$-th root of this. Trying not to connect
   * to terminal, but combine them in one vertex, if possible.
   *
   * The root list is copied from the other SLP.
   *
   * @param other The SLP to compose with *this
   * @return Modified *this
   */
  SLPSet& compose_with(const SLPSet& other);

  //! Check whether or not this program produces the same words as the other one
  /**
   * We check if \p$i\p$-th root in both programs give the same words.
   * Algorithm requires \p$O(n m h)\p$ operations, where \p$n\p$
   * is the size of this->vertices, \p$m\p$ is the size of other.vertices
   * and \p$h\p$ is the height of *this.
   *
   * Note that this function do not reduce the produced words automatically.
   *
   * @param other Program which is compared with *this
   * @return true if programs are equal, false otherwise
   */
  bool equal_to(const SLPSet& other) const;

  //! Reduce words produces by this SLP.
  /**
   * Take each vertex and build a new program such that it produces the freely
   * reduced word.
   *
   * Requires \p$O(n^3 h)\p$ operations, where \p$n\p$ is the number of
   * vertices and  \p$h\p$ is the height of this program.
   *
   * @return The program which produce the freely reduced words.
   */
  SLPSet&& free_reduction() const;
private:

  //! Container for all vertices of the graph.
  std::vector< SLPVertex > vertices;
  std::vector< size_t > roots; //!< Contains the numbers of the vertices in the vertices array

};

#endif	/* COMPOSITIONSYSTEMSET_H */

