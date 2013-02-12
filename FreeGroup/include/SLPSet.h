/*
 * File:   SLPSet.h
 * Author: dpantele
 *
 * Created on November 18, 2012, 4:06 PM
 */

#ifndef SLPSET_H
#define SLPSET_H

#include <gmpxx.h>
typedef mpz_class LongInteger;

#include <vector>
#include <memory>
#include <initializer_list>
#include <iterator>
#include <unordered_map>
#include <functional>
#include <queue>
#include <iostream>

namespace crag {
typedef unsigned int TerminalSymbol;
static const TerminalSymbol INVALID_TERMINAL = 0; //!< Constant representing invalid terminal symbol

//Need this to allow hash access interior of SLPVertex
class SLPVertex;
}

namespace std {
  template<> class hash<crag::SLPVertex>;
}

namespace crag {

//! Internal structure of SLPVertex
struct BasicVertex;

//! One vertex in the straight line program
/**
 * Represents one composition rule of kind \p$A\rightarrow BC \p$.
 * Also stores internal information such as the height of the subtree
 * and the length of it.
 *
 * Should be considered as immutable.
 *
 */
class SLPVertex {
  friend class std::hash<SLPVertex>;
  public:
    //! Default constructor. Constructing 'invalid' vertex.
    SLPVertex()
      : ptr_()
      , negative_(false)
    { }


    static const SLPVertex Null; //Use this vertex to represent invalid vertex

    //Use idiom of "named constructors" to properly create all vertices
    //! Create just terminal vertex.
    static SLPVertex terminal_vertex(const TerminalSymbol&);

    //! Create a new vertex defined by some composition rule
    //TODO: rename it somehow
    static SLPVertex concatenate(const SLPVertex& left_child, const SLPVertex& right_child);

    //Use default copy&move constructors/assignments, not defining them

    //!The implementation of equality operator for #SLPVertex.
    /**
     * Compares the actual addresses of vertices and sign.
     */
    bool operator== (const SLPVertex& other) const {
      return this->ptr_ == other.ptr_ &&
             this->negative_ == other.negative_;
    }

    bool operator!= (const SLPVertex& other) const {
      return !(*this == other);
    }

    //! Returns 'negated' vertex. Null.negate() is Null
    SLPVertex negate() const {
      if (ptr_) {
        return SLPVertex(ptr_, !negative_);
      } else {
        return Null;
      }
    }

    //getters
    inline SLPVertex left_child() const;

    bool has_left_child() const {
      return left_child() != Null;
    }

    inline SLPVertex right_child() const;

    bool has_right_child() const {
      return right_child() != Null;
    }

    inline LongInteger split_point() const {
      return left_child().length();
    }

    inline TerminalSymbol terminal_symbol() const;

    bool is_terminal() const {
      return terminal_symbol() != INVALID_TERMINAL;
    }

    bool is_negative() const {
      return negative_;
    }

    inline LongInteger length() const;

    inline unsigned int height() const;

  protected:
    SLPVertex(const std::shared_ptr<BasicVertex>& ptr, bool negative)
      : ptr_(ptr),
        negative_(negative) {
    }

    std::shared_ptr<BasicVertex> ptr_;  //!< The pointer to SLPVertex
    bool negative_;                     //!< True if we need "inverse" vertex
};



//! Post-order SLP Inspector
/**
 * This is an inspector which goes through the children of some root
 * in the following order: all left children, then all right children,
 * and after that the vertex itself. Useful to produce the resulting
 * words and also to compare SLPs
 */
class SLPPostorderInspector
{
  public:
    //Default constructor is useless, but keep it for #SLPProducedWordIterator()
    SLPPostorderInspector()
    { }

    //! Inspect the subtree of #root
    SLPPostorderInspector(const SLPVertex& root)
      : current_path_({PathState({root, false})}) {
      goto_leftmost_terminal();
    }

    //Use default copy/move assignments/constructors

    //! Get current vertex
    const SLPVertex& current_vertex() const {
      return current_path_.back().vertex;
    }

    //! True if there are no more vertices, we have visited everything
    bool inspection_ended() const {
      return current_path_.empty();
    }

    //! Move inspector to the next vertex
    //TODO: think about this function name
    void go_to_next_vertex();

  private:
    struct PathState {
      SLPVertex vertex;
      bool is_right_child;
    };

    std::vector<PathState> current_path_; //!< Way to the current vertex in the container. We are using this vector like a stack. Each tuple is <Vertex, is_right_child>

    //! Go from the current_path_.back() to the leftmost terminal
    void goto_leftmost_terminal();
};

//! Iterator of #SLPProducedWord.
/**
 * This is a forward iterator for SLPProducedWord. Note that even difference_type is standard ptrdiff_t,
 * it is possible that the actual difference between two iterator may be more than 2^64. However,
 * in all standard algorithms this distance can not be achieved, because it takes to long to move
 * this iterator so long.
 *
 */
class SLPProducedWordIterator : public std::iterator <
        std::forward_iterator_tag,      //iterator_category
        const SLPVertex                 //value_type
  > {
  public:
    SLPProducedWordIterator()
      : inspector_()
      , length_(0)
      , root_()
    { }
    explicit SLPProducedWordIterator(const SLPVertex& root)
      : inspector_(root)
      , length_(0)
      , root_(root)
    { }

    //use default copy/move constructors/assignments

    SLPProducedWordIterator& operator++(); //!< Preincrement
    SLPProducedWordIterator operator++(int) { //!< Postincrement
      SLPProducedWordIterator copy(*this);
      ++(*this);
      return copy;
    }

    reference operator*() const { //!< "Dereference" current symbol
      return inspector_.current_vertex();
    }

    pointer operator->() const {
      return &(inspector_.current_vertex());
    }

    //!< Compare to another iterator.
    /**
     * Compares "#length" && root. If this->length >= root.length() and other->length >= root.length(), then also true
     */
    bool operator==(const SLPProducedWordIterator& other) const {
      return ( root_ == other.root_ &&
               length_ == other.length_
             ) ||
             ( length_ >= root_.length() &&
               other.length_ >= other.root_.length()
             );
    }

    bool operator!=(const SLPProducedWordIterator& other) const {
      return !(*this == other);
    }

  private:
    SLPPostorderInspector inspector_;    //!< Subtree inspector
    LongInteger length_;                 //!< Length already produced
    SLPVertex   root_;                   //!< Root

};

//! Word produced by some #SLPVertex
class SLPProducedWord {
  public:
    typedef SLPVertex value_type;
    typedef SLPProducedWordIterator const_iterator;      //no iterator, only const
    //we do not define size_type, because size() should return LongInteger

    SLPProducedWord() //!< Just empty word
      : root_()
    { }

    explicit SLPProducedWord(const SLPVertex& root) //!< Word produced by some root
      : root_(root)
    { }

   //TODO: define swap, operator ==

    value_type operator[](LongInteger index) const; //!< Get one letter from the word

    const_iterator begin() const { //!< Get the iterator to the first symbol
      return const_iterator(root_);
    }

    const_iterator end() const { //!< Get the iterator to the symbol after the last
      return const_iterator();
    }

    LongInteger size() const {
      return root_.length();
    }

    bool empty() const {
      return size() != 0;
    }

  private:
    SLPVertex root_; //!< The root vertex producing this word
};
}//namespace crag
namespace std {
  //! Definition of the hash for std::pair
  template<typename TFirst, typename TSecond>
  struct hash< std::pair<TFirst, TSecond> > {
  private:
    const std::hash<TFirst> first_hash_;
    const std::hash<TSecond> second_hash_;
  public:
    hash()
      : first_hash_()
      , second_hash_()
    { }
    size_t operator()(const std::pair<TFirst, TSecond>& obj) const {
      size_t first_hash_value = first_hash_(obj.first);
      //Taken from boost/functional/hash
      return second_hash_(obj.second) + 0x9e3779b9 + (first_hash_value << 6) + (first_hash_value >> 2);
    }
  };

  //! Definition of the hash for SignedVertex
  template<>
  struct hash< crag::SLPVertex > {
    private:
      const std::hash<std::shared_ptr<crag::BasicVertex> > ptr_hash_;
    public:
      hash(): ptr_hash_() { }
      size_t operator()(const crag::SLPVertex& vertex) const {
        return vertex.negative_? ~ptr_hash_(vertex.ptr_) : ptr_hash_(vertex.ptr_);
      }
  };
}

namespace crag {
//! Progression tables which are described in the thesis by Lifshits
/**
 * It is a class of the progression table as described in the thesis
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
class SLPMatchingTable {
  public:
    struct MatchResultSequence {
        LongInteger start; //!< The beginning of the first match
        LongInteger step;  //!< The distance between the matches
        LongInteger count; //!< The number of matches

        bool operator==(const MatchResultSequence& other) const {
          return start == other.start && step == other.step && count == other.count;
        }

        bool operator!=(const MatchResultSequence& other) const {
          return !(*this == other);
        }
    };

    const static MatchResultSequence NO_MATCHES;

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
    MatchResultSequence matches(const SLPVertex& pattern,
                                const SLPVertex& text);


  protected:
    std::unordered_map<std::pair<SLPVertex, SLPVertex>, MatchResultSequence> match_table_; //! The actual storage for the calculated values.
};

::std::ostream& operator<<(::std::ostream& os, const SLPMatchingTable::MatchResultSequence& match);

//! class of straight line program collection
/**
 * A straight line program is a program which computes a single word.
 * We used "Polynomial-time Word problems" by Saul Schleimer and 
 * "Algorithms and complexity analysis for processing compressed text"
 * by Yury Lifshits when implementing this class
 *
 * Actually, this class computes several words at the same time, one word for
 * each "root" of this system.
 * TODO make arbitrary number of roots, composition then works under the condition num_terminals_1 == num_roots_2
 */
class SLPSet {
  public:
    SLPSet();
  protected:
    std::vector<SLPVertex> roots;     //!< Contains the roots of this SLP
};

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
 *
 * TODO: Move it to Automorphism class, because it has no sense for abstract SLPS
 */
//SLPSet& compose_with(const SLPSet& other);

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
 *
 * TODO: Put it somewhere, probably into MatchingTable
 */
//bool equal_to(const SLPSet& other) const;

//! Reduce words produces by this SLP.
/**
 * Take each vertex and build a new program such that it produces the freely
 * reduced word.
 *
 * Requires \p$O(n^3 h)\p$ operations, where \p$n\p$ is the number of
 * vertices and  \p$h\p$ is the height of this program. TODO how is it indexed?
 *
 * @return The program which produce the freely reduced words.
 *
 * TODO: Put it into SLPVertex, or maybe some external function, or...
 */
//SLPSet free_reduction() const;

struct BasicVertex {
  public:
    SLPVertex left_child;            //!< Left part of the rule. Use SignedVertex::Null if child is absent.
    SLPVertex right_child;           //!< Right part of the rule.
    TerminalSymbol terminal_symbol;  //!< INVALID_TERMINAL, if non-terminal. Otherwise the number of the symbol, greater that zero.
    LongInteger length;              //!< Length of word produced by the vertex
    unsigned int height;             //!< Height of subtree

    BasicVertex();
};

inline SLPVertex SLPVertex::left_child() const {
  if (!ptr_) {
   return Null;
  }
  if (negative_) {
   return (ptr_->right_child).negate();
  } else {
   return ptr_->left_child;
  }
}

inline SLPVertex SLPVertex::right_child() const {
  if (!ptr_) {
   return Null;
  }
  if (negative_) {
    return ptr_->left_child.negate();
  } else {
    return ptr_->right_child;
  }
}

inline TerminalSymbol SLPVertex::terminal_symbol() const {
  return ptr_ ? ptr_->terminal_symbol : INVALID_TERMINAL;
}

inline LongInteger SLPVertex::length() const {
     return ptr_ ? ptr_->length : 0;
   }

inline unsigned int SLPVertex::height() const {
  return ptr_ ? ptr_->height : 0;
}

//!Namespace for some undocumented internal structures
namespace internal {

//In-order inspector of tree with constraints, required to fill matching table
class SLPMatchingInspector {
  public:
    SLPMatchingInspector() = delete;
    SLPMatchingInspector(const LongInteger& pattern_length, const SLPVertex& text, const LongInteger& left_bound, const LongInteger& interval_length)
      : parent_stack_()
      , pattern_length_(pattern_length)
      , text_position_({text, -left_bound})
      , left_bound_(left_bound)
      , interval_length_(interval_length)
    {
      if (left_bound_ < 0) {
        interval_length_ += left_bound_;
        left_bound_ = 0;
        text_position_.distance_from_left_border = 0;
      }

      if (left_bound_ + interval_length_ > text.length()) {
        interval_length_ = text.length() - left_bound_;
      }

      if (interval_length_ < 0) {
        interval_length_ = 0;
      }
      go_further();
    }

    bool inspection_ended() const {
      return text_position_.vertex == SLPVertex::Null && parent_stack_.empty();
    }

    void go_next() {
      text_position_.distance_from_left_border += text_position_.vertex.left_child().length();
      text_position_.vertex = text_position_.vertex.right_child();
      go_further();
    }

    const SLPVertex& current_text() const {
      return text_position_.vertex;
    }

    LongInteger current_distance_from_left_border() const {
      return text_position_.distance_from_left_border;
    }

    const LongInteger& left_border() const {
      return left_bound_;
    }

    const LongInteger& interval_length() const {
      return interval_length_;
    }

  private:
    struct TextPosition {
      SLPVertex vertex;
      LongInteger distance_from_left_border;
    };
    std::vector<TextPosition> parent_stack_;
    LongInteger pattern_length_;
    TextPosition text_position_;
    LongInteger left_bound_;
    LongInteger interval_length_;

    bool position_is_suitable() const {
      if (text_position_.vertex == SLPVertex::Null ||
          interval_length_ - text_position_.distance_from_left_border < pattern_length_ ||
          text_position_.distance_from_left_border + text_position_.vertex.length()  < pattern_length_
      ) {
        return false;
      }
      return true;
    }

    void go_further();

};

//! Find sequence containing all elements of two sequences and only them
SLPMatchingTable::MatchResultSequence join_sequences(SLPMatchingTable::MatchResultSequence first,
                                                     SLPMatchingTable::MatchResultSequence second);

//! Find sequence containing elements presented in both sequences
SLPMatchingTable::MatchResultSequence intersect_sequences(SLPMatchingTable::MatchResultSequence first,
                                                     SLPMatchingTable::MatchResultSequence second);

//! Find coherent entries of pattern with some bounded inspector
SLPMatchingTable::MatchResultSequence local_search(const SLPVertex& pattern,
                                                   SLPMatchingInspector* inspector,
                                                   SLPMatchingTable* matching_table);
//! SLPMatchingTable::match subroutine
SLPMatchingTable::MatchResultSequence nontrivial_match(const SLPVertex& large_pattern_part,
                                                       const SLPVertex& small_pattern_part,
                                                       bool small_pattern_is_after,
                                                       const SLPVertex& text,
                                                       SLPMatchingTable* matching_table);

} //namespace internal
} //namespace crag
namespace std {
  void swap(crag::SLPMatchingTable::MatchResultSequence& first, crag::SLPMatchingTable::MatchResultSequence& second);
}
#endif  /* SLPSET_H */
