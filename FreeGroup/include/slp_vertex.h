/**
 * \file slp_vertex.h
 * \brief This file defines the basic structure for SLPs. See \ref slp_description for examples.
 */
#pragma once
#ifndef CRAG_FREEGROUP_SLP_VERTEX_H_
#define CRAG_FREEGROUP_SLP_VERTEX_H_

#include <memory>
#include <iostream>
#include <cassert>

#include <gmpxx.h>
typedef mpz_class LongInteger;

#include "boost/pool/pool_alloc.hpp"

namespace crag {
namespace slp {

//!\internal
namespace internal {
class BasicVertex;
struct BasicVertexAllocatorTag{};
typedef boost::fast_pool_allocator<BasicVertex, boost::default_user_allocator_malloc_free, boost::details::pool::null_mutex, 1, 0> BasicVertexPoolAllocator;
}

//! Basic interface for vertex. Also represents empty vertex, use Vertex() as empty vertex
class Vertex {
  public:
    //! Default constructor generating empty vertex.
    constexpr Vertex()
      : vertex_signed_id_(0)
      , vertex_(nullptr)
    {}

    //!Check whether two vertices are equal
    /**
     * Two vertices are considered to be equal if
     *   1. They have the same type
     *   2. They are equal as vertices of that type. See \ref slp_description for details
     *
     * @param[in] other Vertex compared to *this
     * @return true if equal
     */
    bool operator==(const Vertex& other) const {
      return
          vertex_signed_id_ == other.vertex_signed_id_ &&
          //((vertex_ && other.vertex_) || (!vertex_ && !other.vertex_));
          vertex_ == other.vertex_;
    }

    bool operator!=(const Vertex& other) const {
      return !(*this == other);
    }

    //! Return vertex produces the reversed word
    Vertex negate() const {
      return Vertex(-vertex_signed_id_, std::shared_ptr<internal::BasicVertex>(vertex_));
    }

    //! Get the left child of vertex
    inline Vertex left_child() const;

    //! Get the right child of vertex
    inline Vertex right_child() const;

    //! Get the length of left child
    /**
     * Position of the split - place where right child starts. More efficient then
     * left_child().length() since left_child() sometimes need to copy left child.
     *
     * @return The length of left child
     */
    const LongInteger& split_point() const;

    //! Get the length of the word produced by this vertex
    const LongInteger& length() const;

    //! Get the maximum path length to the leaf
    unsigned int height() const;

    //! Internal function to define std::hash<Vertex>
    size_t vertex_hash() const {
      return vertex_id_hash_(vertex_signed_id_);
    }

    //! false only if the vertex is null
    explicit operator bool() const {
      return vertex_signed_id_;
    }

    //! Print debug representation. Used in different tests.
    inline void debug_print(::std::ostream* out) const;

    typedef int64_t VertexSignedId;

    //! This is internal vertex id.
    /**
     * 64-bit integer, sign of it determines whether the vertex is negated,
     * and the absolute value meaning is different:
     *   1. If vertex is terminal, then the value of static_cast<int>(terminal symbol) is stored
     *   2. If vertex is non-terminal, then some auto-incremented value if stored. If the number of nonterminal vertices is greater that \f$2^{63}\f$, then this could case some problems.
     *   3. If vertex is Null, then it is 0.
     *
     * @return Signed vertex id
     */
    VertexSignedId vertex_id() const {
      return vertex_signed_id_;
    }

  protected:
    typedef std::allocator<internal::BasicVertex> VertexAllocator;
    //typedef internal::BasicVertexPoolAllocator VertexAllocator;
    VertexSignedId vertex_signed_id_;
    std::shared_ptr<internal::BasicVertex> vertex_;
    static constexpr std::hash<VertexSignedId> vertex_id_hash_ = std::hash<VertexSignedId>();

    Vertex(VertexSignedId vertex_signed_id, std::shared_ptr<internal::BasicVertex>&& vertex)
      : vertex_signed_id_(vertex_signed_id)
      , vertex_(std::move(vertex))
    {
    }

    static const LongInteger& LongZero();
    static const LongInteger& LongOne();
};

inline void PrintTo(const Vertex& vertex, ::std::ostream* os) {
  vertex.debug_print(os);
}


namespace internal {
class BasicVertex {
  public:
    LongInteger length_;
    unsigned int height_;
    Vertex left_child_;
    Vertex right_child_;

    BasicVertex()
      : length_(0)
      , height_(0)
    { }

    BasicVertex(Vertex&& left_child, Vertex&& right_child)
      : length_(left_child.length() + right_child.length())
      , height_(std::max(left_child.height(), right_child.height()) + 1)
      , left_child_(std::move(left_child))
      , right_child_(std::move(right_child))
    { }
};

}

inline Vertex Vertex::left_child() const {
  if (vertex_) {
    if (vertex_signed_id_ < 0) {
      return vertex_->right_child_.negate();
    } else {
      return vertex_->left_child_;
    }
  } else {
    return Vertex();
  }
}

inline Vertex Vertex::right_child() const {
  if (vertex_) {
    if (vertex_signed_id_ < 0) {
      return vertex_->left_child_.negate();
    } else {
      return vertex_->right_child_;
    }
  } else {
    return Vertex();
  }
}

inline const LongInteger& Vertex::split_point() const {
  if (vertex_) {
    if (vertex_signed_id_ < 0) {
      return vertex_->right_child_.length();
    } else {
      return vertex_->left_child_.length();
    }
  } else {
    return LongZero();
  }
}

inline const LongInteger& Vertex::length() const {
  if (vertex_) {
    return vertex_->length_;
  } else if (vertex_signed_id_) {
    return LongOne();
  } else {
    return LongZero();
  }
}

inline unsigned int Vertex::height() const {
  if (vertex_) {
    return vertex_->height_;
  } else if (vertex_signed_id_) {
    return 1;
  } else {
    return 0;
  }
}

inline void Vertex::debug_print(::std::ostream* out) const {
  if (!vertex_signed_id_) {
    (*out) << "Vertex()";
  } else if (!vertex_) {
    (*out) << "TerminalVertex(" << vertex_signed_id_ << ')';
  } else {
    (*out) << "NonterminalVertex(l=" << vertex_->length_
           << ", h=" << vertex_->height_
           << ", id=" << vertex_signed_id_ << ')';
  }
}


//! Terminal vertex in a SLP. Produces word of length 1.
template <typename TerminalSymbol>
class TerminalVertexTemplate : public Vertex {
  public:
    TerminalVertexTemplate() = delete;

    explicit TerminalVertexTemplate(TerminalSymbol terminal_symbol)
      : Vertex(static_cast<typename Vertex::VertexSignedId>(terminal_symbol), nullptr)
      , terminal_symbol_(terminal_symbol)
    {
      assert(height() == 1 && length() == 1);
    }

    explicit TerminalVertexTemplate(const Vertex& vertex)
      : Vertex(vertex)
      , terminal_symbol_()
    {
      if (vertex_) {
        vertex_ = nullptr;
        vertex_signed_id_ = 0;
      } else {
        terminal_symbol_ = static_cast<TerminalSymbol>(vertex_signed_id_);
      }

      assert((height() == 0 && length() == 0) || (height() == 1 && length() == 1));
    }

    const TerminalSymbol& terminal_symbol() const {
      return terminal_symbol_;
    }

    operator TerminalSymbol() const {
      return terminal_symbol();
    }
  private:
    TerminalSymbol terminal_symbol_;
};

template <typename TerminalSymbol>
::std::ostream& operator << (::std::ostream& stream, const TerminalVertexTemplate<TerminalSymbol>& vertex) {
  return stream << vertex.terminal_symbol();
}

//! Non-terminal vertex in a SLP, represent rule A->BC
class NonterminalVertex : public Vertex {
  public:
    NonterminalVertex() = delete;

    NonterminalVertex(Vertex left, Vertex right)
      : Vertex(
          ++last_vertex_id_ > 0 ? last_vertex_id_ : (last_vertex_id_ = 1),
          std::make_shared<internal::BasicVertex>(
            //get_allocator(),
            ::std::move(left),
            ::std::move(right)
          ))
    {
      assert(left_child());
      assert(right_child());
      assert(left_child().length() > 0);
      assert(right_child().length() > 0);
      assert(height() > 1);
      assert(length() > 1);
    }

  private:
    static Vertex::VertexSignedId last_vertex_id_; //All vertices are enumerated
    static const Vertex::VertexAllocator& get_allocator();
};
}//namespace slp
}//namespace crag

namespace std {
//! Definition of the hash for std::pair
template<typename TFirst, typename TSecond>
struct hash<pair<TFirst, TSecond>> {
private:
  constexpr static hash<TFirst> first_hash_ = hash<TFirst>();
  constexpr static hash<TSecond> second_hash_ = hash<TSecond>();
public:
  size_t operator()(const pair<TFirst, TSecond>& obj) const {
    size_t first_hash_value = first_hash_(obj.first);
    //Taken from boost/functional/hash
    return second_hash_(obj.second) + 0x9e3779b9 + (first_hash_value << 6) + (first_hash_value >> 2);
  }
};

template<typename TFirst, typename TSecond>
constexpr hash<TFirst> hash<pair<TFirst, TSecond>>::first_hash_;
template<typename TFirst, typename TSecond>
constexpr hash<TSecond> hash<pair<TFirst, TSecond>>::second_hash_;


//! Definition of the hash for SignedVertex
template<>
struct hash<crag::slp::Vertex> {
  public:
    size_t operator()(const crag::slp::Vertex& vertex) const {
      return vertex.vertex_hash();
    }
};
} //namespace std


#endif /* CRAG_FREEGROUP_SLP_VERTEX_H_ */
