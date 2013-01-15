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

typedef unsigned int TerminalSymbol;
static const TerminalSymbol INVALID_TERMINAL = 0; //!< Constant representing invalid terminal symbol

#include <vector>
#include <memory>
#include <initializer_list>
#include <iterator>

//! Represents a vertex in program or the "inverse" vertex.
struct SignedVertex {
  public:
    size_t index;     //!< The index of the vertex in in the vertices vector.
    bool negative;    //!< True if we need "inverse" vertex

    SignedVertex()
      : index(-1),
        negative(false)
    { }
};

static const SignedVertex NULL_VERTEX;

bool operator==(const SignedVertex& lhs, const SignedVertex& rhs) {
  return lhs.index == rhs.index &&
      (lhs.negative == rhs.negative ||
       lhs.index == NULL_VERTEX.index);
}

bool operator!=(const SignedVertex& lhs, const SignedVertex& rhs) {
  return !(lhs == rhs);
}

//! Structure representing one vertex in the SLP. Internal.
/**
 * Represents one composition rule of kind \p$A\rightarrow BC \p$.
 * Also stores internal information such as the height of the subtree
 * and the length of it.
 *
 * TODO: check if we can guarantee that any non-terminal vertex has exactly
 * two children.
 *
 */
class SLPVertex {
  public:
    SLPVertex()
        : left_child_(NULL_VERTEX),
          right_child_(NULL_VERTEX),
          terminal_symbol_(INVALID_TERMINAL),
          length_(0),
          height_(0),
          parents_count_(0) {
    }

    /**
     * Get left child of this.
     * @param invert True if inverted
     * @return SignedVertex, left child
     */
    const SignedVertex& left_child(bool invert) const {
      if (invert) {
        return right_child_;
      } else {
        return left_child_;
      }
    }

    bool has_left_child(bool invert) const {
      return left_child(invert) != NULL_VERTEX;
    }

    /**
     * Get right child. @see left_child
     * @param invert True if inverted
     * @return SignedVertex, right child
     */
    const SignedVertex& right_child(bool invert) const {
      return left_child(!invert);
    }

    bool has_right_child(bool invert) const {
      return right_child(invert) != NULL_VERTEX;
    }

    bool is_terminal() const {
      return terminal_symbol_ != INVALID_TERMINAL;
    }

    const TerminalSymbol& terminal_symbol() const {
      return terminal_symbol_;
    }

    const LongInteger& length() const {
      return length_;
    }

    //TODO: add missed getters/setters
  private:
    SignedVertex left_child_; //!< The index of the left vertex in vertices vector. Use CHILD_NOT_EXIST
    SignedVertex right_child_; //!< The index of the right vertex in vertices vector.
    TerminalSymbol terminal_symbol_; //!< NON_TERMINAL, if non-terminal. Otherwise the number of the symbol, greater that zero.
    LongInteger length_;               //!< Length of word produced by the vertex
    unsigned int height_;                   //!< Height of subtree
    unsigned int parents_count_;            //!< Number of parents

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
    MatchResultSequence matches(const SignedVertex& pattern,
                                const SignedVertex& text);

    ProgressionTable(const std::vector<SLPVertex>& pattern_vertices,
                     const std::vector<SLPVertex>& text_vertices)
        : pattern_vertices(pattern_vertices),
          text_vertices(text_vertices),
          table(pattern_vertices.size() * text_vertices.size()) {
    }

  private:
    struct MatchResult {
        MatchResultSequence match; //!< Entries of pattern in text
        MatchResultSequence inversed_match; //!< Entires of inversed pattern in text
    };

    //! Helper function which looks for pattern in text[begin..end]
    std::pair<MatchResultSequence, MatchResultSequence> local_search(
        SignedVertex pattern, SignedVertex text, LongInteger begin, LongInteger end);

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
    SLPSet(unsigned int terminals_count)
        : vertices(2 * terminals_count),
          roots(terminals_count),
          terminals_count(terminals_count) {
      for (size_t terminal_id = 0; terminal_id < terminals_count;
          ++terminal_id) {
        vertices[terminal_id].terminal_symbol_ = terminal_id + 1;
        vertices[terminal_id].parents_count_ = 1;
      }

      for (size_t terminal_id = terminals_count;
          terminal_id < 2 * terminals_count; ++terminal_id) {
        vertices[terminal_id].length_ = 1;
        vertices[terminal_id].height_ = 1;
        vertices[terminal_id].left_child_ = SignedVertex(
            terminal_id - terminals_count, false);
        roots[terminal_id - terminals_count] = terminal_id;
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
     * @param rule_rhs The right side of prudction rule, see docs for rule_lhs.
     *                 if rule_rhs == 0, then the second child is epsont.
     */

    SLPSet(unsigned int terminals_count, unsigned int nontrivial_root,
           int rule_lhs, int rule_rhs)
        : SLPSet(terminals_count) {
      SLPVertex non_trivial;

      if (rule_lhs != -rule_rhs) { //If one vertex cancels another one, leave just root without any children
        non_trivial.left_child_.index = abs(rule_lhs) - 1;
        non_trivial.left_child_.is_negative = (rule_lhs < 0);
        non_trivial.height_ = 1;

        if (rule_rhs != 0) {
          non_trivial.right_child_.index = abs(rule_rhs) - 1;
          non_trivial.right_child_.is_negative = (rule_rhs < 0);
          non_trivial.length_ = 2;
        } else {
          non_trivial.length_ = 1;
        }

        if (abs(rule_lhs) != nontrivial_root) {
          vertices[non_trivial.left_child_.index].parents_count_ += 1;
        }

        if (abs(rule_rhs) != nontrivial_root) {
          vertices[non_trivial.right_child_.index].parents_count_ += 1;
        }
      }

      vertices[nontrivial_root + terminals_count - 1] = non_trivial;
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

    //! Post-order SLP Inspector
    /**
     * This is an inspector which goes through the
     * vertices of the SLP in the following orders: first, all left children,
     * then all right children, and after that the vertex itself. Useful to
     * produce the resulting words and also to compare SLPs
     */
    class ConstVerticesPostorderInspector
    {
      public:
      ConstVerticesPostorderInspector()
      : vertices(nullptr)
      {}
      //! Inspect the subtree of #root in #vertices
      ConstVerticesPostorderInspector(const std::vector<SLPVertex>& vertices, const SignedVertex& root)
      : vertices(&vertices)
      {
        this->current_path.push(root);
        this->goto_leftmost_terminal();
      }

      explicit ConstVerticesPostorderInspector(const ConstVerticesPostorderInspector& other)
      : vertices(other.vertices)
      , current_path(other.current_path)
      {}

      ConstVerticesPostorderInspector& operator=(const ConstVerticesPostorderInspector& other) {
        this->vertices = other.vertices;
        this->current_path = other.current_path;

        return *this;
      }

      //! Return the index of current vertex
      const SignedVertex& current_vertex() const {
        return current_path.top();
      }

      //! True if there are no more vertices, we have visited everything
      bool inspection_ended() const {
        return this->current_path.empty();
      }

      //! Move inspector to the next vertex
      void next_vertex() {
        SignedVertex current = this->current_path.top();
        this->current_path.pop();

        if(this->inspection_ended()) {
          return;
        }
        SignedVertex parent = this->current_path.top();
        SignedVertex right_sibling = SLPSet::get_right_child(parent, *(this->vertices));

        if (right_sibling == current || right_sibling == SLPVertex::CHILD_NOT_EXIST) {
          //We have already visited all right children of the parent vertex, stop on parent
        } else {
          //We have visited all left children of parent, going to the right
          this->current_path.push(right_sibling);
          this->goto_leftmost_terminal();
        }
      }

      private:
      const std::vector< SLPVertex >* vertices; //!< Pointer to vertices container
      std::stack< SignedVertex > current_path;//!< Way to the current vertex in the container

      //!Go from the current_path->top() to the leftmost terminal
      void goto_leftmost_terminal() {
        while (this->current_path.top() != SLPVertex::CHILD_NOT_EXIST) {
          bool inversed = this->current_path.top().is_negative;
          const SignedVertex& left_child = SLPSet::get_left_child(this->current_path.top(), *(this->vertices));
          if (left_child != SLPVertex::CHILD_NOT_EXIST) {
            this->current_path.push(left_child);
          } else {
            //Don't have any left children, then go right
            this->current_path.push(SLPSet::get_right_child(this->current_path.top(), *(this->vertices)));
          }

          this->current_path.top().is_negative ^= inversed;
        }

        this->current_path.pop();
      }
    };

    private:

    //! Container for all vertices of the graph.
    std::vector<SLPVertex> vertices;
    std::vector<size_t> roots;//!< Contains the numbers of the vertices in the vertices array
    unsigned int terminals_count;//!< The number of terminals in the system.

    //! Helper function to get left child with respect to the sign of vertex
    static const SignedVertex& get_left_child(const SignedVertex& vertex, const std::vector< SLPVertex >& vertices) {
      if (vertex.is_negative) {
        return vertices[vertex.index].right_child_;
      } else {
        return vertices[vertex.index].left_child_;
      }
    }

    //! Helper function to get right child with respect to the sign of vertex
    static const SignedVertex& get_right_child(const SignedVertex& vertex, const std::vector< SLPVertex >& vertices) {
      if (! vertex.is_negative) {
        return vertices[vertex.index].right_child_;
      } else {
        return vertices[vertex.index].left_child_;
      }
    }

  };

//!The standard-interfaced iterator for SLPProducedWord class
class SLPProducedWordConstIterator: public std::iterator<
    std::forward_iterator_tag,      //iterator_category
    const TerminalSymbol,           //value_type
    LongInteger                     //difference_type
> {
  public:
    SLPProducedWordConstIterator()
        : inspector(), vertices(nullptr), current_word_length(0) {
    }

    SLPProducedWordConstIterator(const std::vector<SLPVertex>& vertices,
                                 const SignedVertex& root)
        : inspector(vertices, root), vertices(&vertices), current_word_length(0) {
    }

    explicit SLPProducedWordConstIterator(
        const SLPProducedWordConstIterator& other)
        : inspector(other.inspector),
          current_word_length(other.current_word_length) {
    }

    SLPProducedWordConstIterator& operator=(
        const SLPProducedWordConstIterator& other) {
      inspector = other.inspector;
      current_word_length = other.current_word_length;
    }

    bool operator==(const SLPProducedWordConstIterator& other) {
      return vertices == (other.vertices)
          && (vertices == nullptr || current_word_length
              == other.current_word_length);
    }

    bool operator!=(const SLPProducedWordConstIterator& other) {
      return !(*this == other);
    }

    SLPProducedWordConstIterator& operator++() {
      inspector.next_vertex();

      while (vertices->at(inspector.current_vertex().index).terminal_symbol_
          == SLPVertex::NON_TERMINAL) {
        inspector.next_vertex();
      }

      ++current_word_length;
    }

    const TerminalSymbol& operator*() {
      return (inspector.current_vertex().is_negative ? -1 : 1)
          * vertices->at(inspector.current_vertex().index).terminal_symbol_;
    }

  private:
    SLPSet::ConstVerticesPostorderInspector inspector; //!<We use this inspector to list all terminals
    const std::vector<SLPVertex>* vertices;         //!< Container with vertices
    LongInteger current_word_length; //!<The current position in the resulting word
};

//! Class which represents a constant word produced by some vertex from SLP
/**
 * This class supports inspections of
 *
 */
class SLPProducedConstWord {
  public:
    //! The only constructor for this word
    /**
     * Construct a word adapter on SLP specified by vertices vector with root vertex
     * having index root_vertex in the vector of vertices.
     *
     * @param vertices The vector of vertices in SLP
     * @param root_vertex The index of the root vertex in SLP
     */
    SLPProducedConstWord(const std::vector<SLPVertex>& vertices,
                         size_t root_vertex)
        : vertices(vertices), root_vertex(root_vertex) {
    }

    unsigned int get_letter(size_t index) const;
    size_t size() const;

  private:
    SLPProducedConstWord() {
    } //word without associated vertices is useless

    const std::vector<SLPVertex>& vertices;
    size_t root_vertex;
};

#endif	/* COMPOSITIONSYSTEMSET_H */

