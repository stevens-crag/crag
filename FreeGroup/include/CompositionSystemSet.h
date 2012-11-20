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

//! class of composition system collection
/**
 * A composition system is a program which computes a single word.
 * We used "Polynomial-time Word problems" by Saul Schleimer when implementing
  this class.
 *
 * Actually, this class computes several words at the same time, one word for
 * each "root" of this system.
 */
class CompositionSystemSet {
public:
  //! Default constructor
  CompositionSystemSet();

  //! Default copy constructor
  explicit CompositionSystemSet(const CompositionSystemSet& other);

  //! Default virtual destructor
  virtual ~CompositionSystemSet() {}

  //! Compose current system with another one
  /**
   * Add all vertices from the other system to *this and connect \p$i\p$-th
   * terminal of other system to \p$i\p$-th root of this. Trying not to connect
   * to terminal, but combine them in one vertex, if possible.
   *
   * The root list is copied from the other system.
   *
   * @param other The composition system to combine with *this
   * @return Modified *this
   */
  CompositionSystemSet& composeWith(const CompositionSystemSet& other);
private:

  //! Struct representing one vertex in the compositions system. Internal.
  /**
   * Represents one composition rule of kind \p$A\rightarrow B[i:j]C[k:l] \p$.
   * Also stores internal information such as the height of the subtree
   * and the length of it.
   *
   * TODO: check if we can gurantee that any non-terminal vertex has exactly
   * two children.
   *
   * TODO: maybe store the
   */
  struct Vertex {
    int left_child;                        //!< The index of the left vertex in vertices vector.
    LongInteger left_child_left_boundary;  //!< The number of letters to throw from the beginning of the result. Use 0 to keep everything.
    LongInteger left_child_right_boundary; //!< The result will end just before getting to this boundary. If it is less that the left_boundary, then no limit is set. Use -1 to keep everything.
    int right_child;                       //!< The index of the right vertex in vertices vector.
    LongInteger right_child_left_boundary; 
    LongInteger right_child_right_boundary;
    bool inverted;                         //!< Produce \p$ C^{-1} B^{-1} \p$ instead od \p$ B C\p$

    //! Get the length of the word produced by subtree with this vertex as a root.
    inline LongInteger getLength() {
      return left_child_right_boundary - left_child_left_boundary +
              right_child_right_boundary - right_child_left_boundary;
    }

    Vertex()
      : left_child(-1)
      , left_child_left_boundary(-1)
      , left_child_right_boundary(-1)
      , right_child(-1)
      , right_child_left_boundary(-1)
      , right_child_right_boundary(-1)
      , inverted(false)
    { }
  };

  //! Container for all vertices of the graph.
  /**
   * This container is shared between different instances of CompostitionSystemSet
   * as much as possible by storing different root for each collection. 
   */
  std::shared_ptr< std::vector< Vertex > > vertices;

  std::vector<size_t> roots; //!< Contains the numbers of the vertices in the vertice.s array

};


#endif	/* COMPOSITIONSYSTEMSET_H */

