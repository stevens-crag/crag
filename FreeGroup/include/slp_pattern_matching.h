/*
 * slp_pattern_matching.h
 *
 *  Created on: Feb 26, 2013
 *      Author: dpantele
 */

#pragma once
#ifndef CRAG_FREEGROUP_SLP_PATTERN_MATCHING_H_
#define CRAG_FREEGROUP_SLP_PATTERN_MATCHING_H_

#include <unordered_map>
#include <iostream>

#include "slp_vertex.h"
#include "slp_inspector.h"
#include "arithmetic_sequence.h"

namespace crag {
namespace slp {

//! Matching tables (or progression tables) which are described in the thesis by Lifshits
/**
 * It is a class of the matching table as described in the thesis
 * by Yury Lifshits. This table enumerates the entries of the subtrees
 * of pattern to the subtrees of text.
 *
 * It is calculated for two SLPs \em P and \em T, and has \em nm cells, where
 * \em n and \em are the numbers of vertices in \em P and  \em T correspondingly.
 * In the cell PT[i][j] we have all entries of the word produced by
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
 */
class MatchingTable {
  public:
    //! Return all matches of the pattern around the "split point"
    /**
     * This function get the result from #match_table_ and recursively calculate
     * it if needed.
     *
     * @param pattern The SLP for a word we want to find in text.
     * @param text    The SLP for the text where we are searching.
     *
     * @return The sequence of the beginnings of matches.
     */
    const FiniteAritmeticSequence& matches(const Vertex& pattern,
                                           const Vertex& text);

  protected:
    ::std::unordered_map<std::pair<Vertex, Vertex>, FiniteAritmeticSequence> match_table_; //! The actual storage for the calculated values.
};

namespace inspector {

class BoundedTaskAcceptor : public TaskAcceptor {
  public:
    const LongInteger& first_lookup_end_position_;
    const LongInteger& last_lookup_begin_position_;

    BoundedTaskAcceptor(const LongInteger& first_lookup_end_position, const LongInteger& last_lookup_begin_position)
      : TaskAcceptor()
      , first_lookup_end_position_(first_lookup_end_position)
      , last_lookup_begin_position_(last_lookup_begin_position)
    { }

    bool accept(const InspectorTask& task) {
       return TaskAcceptor::accept(task) &&
           task.left_siblings_length + task.vertex.length() >= first_lookup_end_position_ && //current vertex should fit pattern
           task.left_siblings_length <= last_lookup_begin_position_; //
    }
};

};

class PatternMatchesGenerator {
  public:
    PatternMatchesGenerator(const Vertex& pattern, const Vertex& text, LongInteger lookup_from, LongInteger lookup_length, ::std::shared_ptr<MatchingTable> matching_table)
      : first_lookup_begin_position_(lookup_from)
      , first_lookup_end_position_(lookup_from + pattern.length())
      , last_lookup_begin_position_(((lookup_from += lookup_length) > text.length()) ? text.length() - pattern.length() : lookup_from - pattern.length())
      , text_inspector_(text, ::std::unique_ptr<inspector::BoundedTaskAcceptor>(new inspector::BoundedTaskAcceptor(
        first_lookup_end_position_, last_lookup_begin_position_
      )))
      , pattern_(pattern)
      , text_(text)
      , matching_table_(matching_table)
    {
      std::cout << "first_lookup_begin_position_: " << first_lookup_begin_position_ << std::endl;
      std::cout << "first_lookup_end_position_: " << first_lookup_end_position_ << std::endl;
      std::cout << "last_lookup_begin_position_: " << last_lookup_begin_position_ << std::endl;
      std::cout << "4: " << lookup_from << std::endl;
      std::cout << "5: " << lookup_length << std::endl;
      std::cout << "6: " << pattern.length() << std::endl;
      std::cout << "7: " << text.length() << std::endl;

    }

    PatternMatchesGenerator(const Vertex& pattern, const Vertex& text, LongInteger lookup_from, LongInteger lookup_length)
      : PatternMatchesGenerator(pattern, text, lookup_from, lookup_length, ::std::make_shared<MatchingTable>())
    {
      next_match();
    }

    const FiniteAritmeticSequence& next_match();

    const FiniteAritmeticSequence& current_match() const {
      return current_match_;
    }
  protected:
    LongInteger first_lookup_begin_position_;
    LongInteger first_lookup_end_position_;
    LongInteger last_lookup_begin_position_;
    InorderInspector text_inspector_;
    const Vertex& pattern_;
    const Vertex& text_;
    FiniteAritmeticSequence current_match_;
    ::std::shared_ptr<MatchingTable> matching_table_;
};


} //namespace slp
} //namespace crag


#endif /* CRAG_FREEGROUP_SLP_PATTERN_MATCHING_H_ */
