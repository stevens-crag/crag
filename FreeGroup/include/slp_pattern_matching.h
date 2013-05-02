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
    const FiniteArithmeticSequence& matches(const Vertex& pattern,
                                           const Vertex& text);

    bool is_calculated(const Vertex& pattern, const Vertex& text) const {
      return pattern.length() <= 1 || text.length() <= 1 || match_table_->count(::std::make_pair(pattern, text));
    }

    size_t size() const {
      return match_table_->size();
    }

    MatchingTable()
      : match_table_(::std::make_shared<std::unordered_map<std::pair<Vertex, Vertex>, FiniteArithmeticSequence>>())
    { }

    MatchingTable clone() const {
      return MatchingTable(::std::make_shared<std::unordered_map<std::pair<Vertex, Vertex>, FiniteArithmeticSequence>>(*match_table_));
    }

  protected:
    MatchingTable(::std::shared_ptr<std::unordered_map<std::pair<Vertex, Vertex>, FiniteArithmeticSequence>> match_table)
      : match_table_(::std::move(match_table))
    { }
    ::std::shared_ptr<std::unordered_map<std::pair<Vertex, Vertex>, FiniteArithmeticSequence>> match_table_; //! The actual storage for the calculated values.
};

namespace inspector {

class BoundedTaskAcceptor {
  public:
    const LongInteger& first_lookup_end_position_;
    const LongInteger& last_lookup_begin_position_;
    const LongInteger& pattern_length_;

    BoundedTaskAcceptor(const LongInteger& first_lookup_end_position, const LongInteger& last_lookup_begin_position, const LongInteger& pattern_length)
      : first_lookup_end_position_(first_lookup_end_position)
      , last_lookup_begin_position_(last_lookup_begin_position)
      , pattern_length_(pattern_length)
    { }

    bool operator()(const InspectorTask& task) const {
       return task.vertex.length() >= pattern_length_ &&
           task.left_siblings_length + task.vertex.length() >= first_lookup_end_position_ &&
           task.left_siblings_length <= last_lookup_begin_position_;
    }
};

};

class PatternMatchesGenerator {
  public:
    PatternMatchesGenerator(const Vertex& pattern, const Vertex& text, LongInteger lookup_from, LongInteger lookup_length, MatchingTable* matching_table)
      : first_lookup_begin_position_(lookup_from)
      , first_lookup_end_position_(lookup_from + pattern.length())
      , last_lookup_begin_position_(((lookup_from += lookup_length) > text.length()) ? text.length() - pattern.length() : lookup_from - pattern.length())
      , text_inspector_(text, inspector::BoundedTaskAcceptor(
          first_lookup_end_position_, last_lookup_begin_position_, pattern.length()
        ))
      , pattern_(pattern)
      , text_(text)
      , matching_table_(matching_table ? *matching_table : MatchingTable())
    { }


    PatternMatchesGenerator(const Vertex& pattern, const Vertex& text, LongInteger lookup_from, LongInteger lookup_length)
      : PatternMatchesGenerator(
          pattern,
          text,
          ::std::move(lookup_from),
          ::std::move(lookup_length),
          nullptr
        )
    { }

    FiniteArithmeticSequence next_match();

    explicit operator bool() const {
      return !text_inspector_.stopped();
    }
  protected:
    LongInteger first_lookup_begin_position_;
    LongInteger first_lookup_end_position_;
    LongInteger last_lookup_begin_position_;
    Inspector<inspector::Inorder, inspector::BoundedTaskAcceptor> text_inspector_;
    const Vertex& pattern_;
    const Vertex& text_;
    MatchingTable matching_table_;
};

namespace internal {
FiniteArithmeticSequence nontrivial_match(
    const Vertex& large_pattern_part,
    const Vertex& small_pattern_part,
    bool small_pattern_is_after,
    const Vertex& text,
    MatchingTable* matching_table);

} //namespace internal

} //namespace slp
} //namespace crag


#endif /* CRAG_FREEGROUP_SLP_PATTERN_MATCHING_H_ */
