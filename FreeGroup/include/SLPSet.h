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

//! Struct representing one vertex in the SLP. Internal.
/**
 * Represents one composition rule of kind \p$A\rightarrow BC \p$.
 * Also stores internal information such as the height of the subtree
 * and the length of it.
 *
 * TODO: check if we can gurantee that any non-terminal vertex has exactly
 * two children.
 *
 */
struct SLPVertex {
  const size_t CHILD_NOT_EXIST = -1;     //!< Constant indicating that there is no such child.
  size_t left_child;                     //!< The index of the left vertex in vertices vector. Use CHILD_NOT_EXIST
  size_t right_child;                    //!< The index of the right vertex in vertices vector.
  bool inverted;                         //!< Produce \p$ C^{-1} B^{-1} \p$ instead of \p$ B C\p$
  int terminal_symbol;                   //!< NON_TERMINAL, if non-terminal. Otherwise the number of the symbol.
  const int NON_TERMINAL = 0;            //!< Constant, meaning this vertex is non-terminal
  LongInteger length;                    //!< Length of word produced by the vertex
  unsigned int parents_counter;          //!< Number of parents

  CompositionSystemVertex()
    : left_child(CHILD_NOT_EXIST)
    , right_child(CHILD_NOT_EXIST)
    , inverted(false)
    , terminal_symbol(NON_TERMINAL)
    , length(0)
  { }
};

//! class of composition system collection
/**
 * A composition system is a program which computes a single word.
 * We used "Polynomial-time Word problems" by Saul Schleimer when implementing
  this class.
 *
 * Actually, this class computes several words at the same time, one word for
 * each "root" of this system.
 */
class SLPSet {
public:
  //! Default constructor
  SLPSet();

  //! Default copy constructor
  explicit SLPSet(const SLPSet& other);

  //! Default virtual destructor
  virtual ~SLPSet() {}

 //! Compose current system with another one
  /**
   * Add all vertices from the other system to *this and connect \p$i\p$-th
   * terminal of other system to \p$i\p$-th root of this. Trying not to connect
   * to terminal, but combine them in one vertex, if possible.
   *
   * The root list is copied from the other system.
   *
   * @param other The composition system to compose with *this
   * @return Modified *this
   */
  SLPSet& compose_with(const SLPSet& other);

  bool equal_to(const SLPSet& other) const;

  SLPSet free_reduction() const;
private:

  //! Container for all vertices of the graph.
  std::vector< SLPVertex > vertices;
  std::vector<size_t> roots; //!< Contains the numbers of the vertices in the vertices array

  struct SignedVertex {
    size_t index;
    bool is_negative;

    SignedVertex(size_t index, bool negative)
            : index(index)
            , is_negative(negative)
    { }
  };

  struct PatternMatch {
    LongInteger start;
    LongInteger step;
    LongInteger count;
  };

  void compute_PT(SignedVertex pattern,
                  SignedVertex text,
                  std::vector< std::vector< PatternMatch > > * prefix_table
  );

  void local_search(SignedVertex pattern, SignedVertex text, LongInteger begin, LongInteger end);
};
#endif	/* COMPOSITIONSYSTEMSET_H */

