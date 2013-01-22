/*
 * File:   SLPSet.h
 * Author: dpantele
 *
 * Created on November 18, 2012, 4:06 PM
 */

#ifndef SLPSET_H
#define	SLPSET_H

#include <gmpxx.h>
typedef mpz_class LongInteger;

#include <vector>
#include <memory>
#include <initializer_list>
#include <iterator>
#include <unordered_map>
#include <functional>
#include <queue>

using std::hash;

typedef unsigned int TerminalSymbol;
static const TerminalSymbol INVALID_TERMINAL = 0; //!< Constant representing invalid terminal symbol

//Need this to allow hash access interior of SLPVertex
class SLPVertex;
namespace std {
  template<> class hash<SLPVertex>;
}

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

    inline TerminalSymbol terminal_symbol() const;

    bool is_terminal() const {
      return terminal_symbol() != INVALID_TERMINAL;
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
      : current_path_() {
      current_path_.push_back(root);
      goto_leftmost_terminal();
    }

    //Use default copy/move assignments/constructors

    //! Get current vertex
    const SLPVertex& current_vertex() const {
      return current_path_.back();
    }

    //! True if there are no more vertices, we have visited everything
    bool inspection_ended() const {
      return current_path_.empty();
    }

    //! Move inspector to the next vertex
    //TODO: think about this function name
    void go_to_next_vertex();

  private:
    std::vector< SLPVertex > current_path_; //!< Way to the current vertex in the container. We are using this vector like a stack

    //! Go from the current_path_.back() to the leftmost terminal
    void goto_leftmost_terminal();
};

/* We can't emulate standard container here to use functions from STL, because
 * the difference between iterators must be integral type, which is impossible here.
 *
 * So it will be -like interfaces.
 */
//! Iterator-like interface of #SLPProducedWord.
/**
 * Due to the fact that difference between iterators is the LongIntger,
 * we decided not to implement STL iterator interface, so it is just iterator-like
 * implementation.
 *
 */
class SLPProducedWordIterator {
  public:
    SLPProducedWordIterator()
      : inspector_()
      , length_(0)
      , root_()
    { }
    explicit SLPProducedWordIterator(const SLPVertex& root);

    SLPProducedWordIterator& operator++(); //!< Prefix increment
    SLPProducedWordIterator operator++(int); //!< Postfix increment

    const TerminalSymbol& operator*() const; //!< "Dereference" current symbol

    //!< Compare to another iterator.
    /**
     * Compares "#length" && root. If this->length >= root.length() and other->length >= root.length(), then also true
     */
    bool operator==(const SLPProducedWordIterator& other);
    bool operator!=(const SLPProducedWordIterator& other);

  private:
    SLPPostorderInspector inspector_; //!< Subtree inspector
    LongInteger length_;              //!< Length already produced
    SLPVertex   root_;                //!< Root

};

//! Word produced by some #SLPVertex
class SLPProducedWord {
  public:
    SLPProducedWord(); //!< Just empty word
    explicit SLPProducedWord(const SLPVertex& root); //!< Word produced by some root

    const TerminalSymbol& operator[](LongInteger index) const; //!< Get one letter from the word

    SLPProducedWordIterator begin() const; //!< Get the iterator to the first symbol
    SLPProducedWordIterator end() const;   //!< Get the iterator to the symbol after the last
  private:
    SLPVertex root_; //!< The root vertex producing this word
};

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
  struct hash< SLPVertex > {
    private:
      const std::hash<std::shared_ptr<BasicVertex> > ptr_hash_;
    public:
      hash(): ptr_hash_() { }
      size_t operator()(const SLPVertex& vertex) const {
        return vertex.negative_? ~ptr_hash_(vertex.ptr_) : ptr_hash_(vertex.ptr_);
      }
  };
}

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
    MatchResultSequence matches(const SLPVertex& pattern,
                                const SLPVertex& text);


  protected:
    //! Helper function which looks for pattern in text[begin..end]
    std::pair<MatchResultSequence, MatchResultSequence> local_search(
        const SLPVertex&  pattern, const SLPVertex&  text, LongInteger begin, LongInteger end);

    std::unordered_map<std::pair<SLPVertex, SLPVertex>, MatchResultSequence> match_table_; //! The actual storage for the calculated values.
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

#endif	/* SLPSET_H */

