/*
 * slp_vertex.h
 *
 *  Created on: Feb 10, 2013
 *      Author: dpantele
 */

#ifndef CRAG_FREEGROUP_SLP_VERTEX_H_
#define CRAG_FREEGROUP_SLP_VERTEX_H_

#include <gmpxx.h>
typedef mpz_class LongInteger;

#include <memory>

namespace crag {
namespace slp {

class Vertex;

namespace internal {
class BasicVertex {
  public:
    LongInteger length_;
    unsigned int height_;

    BasicVertex()
      : length_(0)
      , height_(0)
    { }

    BasicVertex(LongInteger&& length, unsigned int height)
      : length_(::std::move(length))
      , height_(height)
    { }

    virtual ~BasicVertex() {}

    virtual Vertex left_child() const;
    virtual Vertex right_child() const;
    virtual Vertex negate() const;

    virtual size_t vertex_hash() const {
      return 0;
    }

    virtual bool call_is_same_vertex(const BasicVertex& other) const {
      return true;
    }

    bool is_same_vertex(const BasicVertex& other) const {
      return true;
    }
};
}

//! Basic interface for vertex. Also represents empty vertex, use Vertex::Null as empty vertex
class Vertex {
  public:
    const static Vertex Null; //!< Use this vertex to represent invalid (or empty) vertex

    bool operator==(const Vertex& other) const {
      return (!vertex_) ?
          !other.vertex_ : other.vertex_ &&
          typeid(*vertex_) == typeid(*other.vertex_) && vertex_->call_is_same_vertex(*other.vertex_);
    }

    bool operator!=(const Vertex& other) const {
      return !(*this == other);
    }

    //! Return vertex produces the reversed word
    Vertex negate() const {
      return vertex_ ? vertex_->negate() : Null;
    }

    Vertex left_child() const {
      return vertex_ ? vertex_->left_child() : Null;
    }

    Vertex right_child() const {
      return vertex_ ? vertex_->right_child() : Null;
    }

    const LongInteger& length() const {
      return vertex_ ? vertex_->length_ : LongZero;
    }

    unsigned int height() const {
      return vertex_ ? vertex_->height_ : 0;
    }

    size_t vertex_hash() const {
      return vertex_ ? vertex_->vertex_hash() : 0;
    }

  protected:
    ::std::shared_ptr<internal::BasicVertex> vertex_;
    Vertex() {} //!< Default constructor generating empty vertex. Use #Null instead of it
    Vertex(::std::shared_ptr<internal::BasicVertex>&& vertex)
      : vertex_(vertex)
    { }
    const static LongInteger LongZero; //We use it to return length of Null vertex

};

const Vertex Vertex::Null;
const LongInteger Vertex::LongZero;

inline Vertex internal::BasicVertex::negate() const {
  return Vertex::Null;
}

inline Vertex internal::BasicVertex::left_child() const {
  return Vertex::Null;
}

inline Vertex internal::BasicVertex::right_child() const {
  return Vertex::Null;
}


namespace internal {

template <typename TerminalSymbol>
class BasicTerminalVertex : public BasicVertex {
  public:
    TerminalSymbol terminal_symbol_;

    virtual Vertex negate() const;

    virtual size_t vertex_hash() const {
      return terminal_symbol_hash(terminal_symbol_);
    }

    virtual bool call_is_same_vertex(const BasicVertex& other) const {
      return is_same_vertex(dynamic_cast<const BasicTerminalVertex&>(other));
    }

    bool is_same_vertex(const BasicTerminalVertex& other) const {
      return terminal_symbol_ == other.terminal_symbol_;
    }

    BasicTerminalVertex(const TerminalSymbol& terminal_symbol)
      : BasicVertex(LongInteger(1), 1)
      , terminal_symbol_(terminal_symbol)
    { }

  private:
    const static ::std::hash<TerminalSymbol> terminal_symbol_hash;
};
}


//! Terminal vertex in a SLP. Produces word of length 1.
template <typename TerminalSymbol>
class TerminalVertexTemplate : public Vertex {
  public:
    TerminalVertexTemplate() = delete;

    explicit TerminalVertexTemplate(const TerminalSymbol& terminal_symbol)
      : Vertex(::std::make_shared<internal::BasicTerminalVertex<TerminalSymbol>>(terminal_symbol))
      , terminal_vertex_ptr_(static_cast<internal::BasicTerminalVertex<TerminalSymbol>*>(vertex_.get()))
    { }

    explicit TerminalVertexTemplate(const Vertex& vertex)
      : Vertex(vertex)
      , terminal_vertex_ptr_(dynamic_cast<internal::BasicTerminalVertex<TerminalSymbol>*>(vertex_.get()))
    {
      if (!terminal_vertex_ptr_) {
        vertex_ = 0;
      }
    }

    const TerminalSymbol& terminal_symbol() const {
      return terminal_vertex_ptr_ ? terminal_vertex_ptr_->terminal_symbol_ : NullTerminalSymbol;
    }

    operator TerminalSymbol() const {
      return terminal_symbol();
    }

  private:
    const internal::BasicTerminalVertex<TerminalSymbol>* terminal_vertex_ptr_;
    const static TerminalSymbol NullTerminalSymbol;

};

template <typename TerminalSymbol>
::std::ostream& operator << (::std::ostream& stream, const TerminalVertexTemplate<TerminalSymbol>& vertex) {
  return stream << vertex.terminal_symbol();
}

template <typename TerminalSymbol>
Vertex internal::BasicTerminalVertex<TerminalSymbol>::negate() const {
  return TerminalVertexTemplate<TerminalSymbol>(-terminal_symbol_);
}

template <typename TerminalSymbol>
const ::std::hash<TerminalSymbol> internal::BasicTerminalVertex<TerminalSymbol>::terminal_symbol_hash = ::std::hash<TerminalSymbol>();

template <typename TerminalSymbol>
const TerminalSymbol TerminalVertexTemplate<TerminalSymbol>::NullTerminalSymbol = TerminalSymbol();

namespace internal {
struct NonterminalVertexNodeData {
  Vertex left_child;
  Vertex right_child;
};

class BasicNonterminalVertex : public internal::BasicVertex {
  public:
    ::std::shared_ptr<NonterminalVertexNodeData> node_data_ptr_;
    bool negate_node_;
    static const ::std::hash<std::shared_ptr<NonterminalVertexNodeData>> ptr_hash;

    virtual size_t vertex_hash() const {
      return negate_node_ ? ~ptr_hash(node_data_ptr_) : ptr_hash(node_data_ptr_);
    }

    virtual Vertex negate() const;

    virtual Vertex left_child() const {
      return negate_node_ ? node_data_ptr_->right_child.negate() : node_data_ptr_->left_child;
    }

    virtual Vertex right_child() const {
      return negate_node_ ? node_data_ptr_->left_child.negate() : node_data_ptr_->right_child;
    }

    virtual bool call_is_same_vertex(const BasicVertex& other) const {
      return is_same_vertex(dynamic_cast<const BasicNonterminalVertex&>(other));
    }

    bool is_same_vertex(const BasicNonterminalVertex& other) const {
      return node_data_ptr_ == other.node_data_ptr_ && negate_node_ == other.negate_node_;
    }

    BasicNonterminalVertex(::std::shared_ptr<NonterminalVertexNodeData>&& node_data_ptr, bool negate_node)
      : BasicVertex(
          node_data_ptr->left_child.length() + node_data_ptr->right_child.length(),
          ::std::max(node_data_ptr->left_child.height(), node_data_ptr->left_child.height()) + 1
        )
      , node_data_ptr_(std::move(node_data_ptr))
      , negate_node_(negate_node)
    { }
};
}

const ::std::hash<std::shared_ptr<internal::NonterminalVertexNodeData>> internal::BasicNonterminalVertex::ptr_hash;

//! Non-terminal vertex in a SLP, represent rule A->BC
class NonterminalVertex : public Vertex {
  public:
    NonterminalVertex() = delete;

    template <typename LeftVertexT, typename RightVertexT>
    NonterminalVertex(LeftVertexT&& left, RightVertexT&& right)
      : Vertex(::std::make_shared<internal::BasicNonterminalVertex>(
          ::std::make_shared<internal::NonterminalVertexNodeData>(internal::NonterminalVertexNodeData({left, right})),
          false
        ))
    { }

  private:
    friend class internal::BasicNonterminalVertex;

    NonterminalVertex(::std::shared_ptr<internal::BasicNonterminalVertex>&& vertex)
      : Vertex(std::move(vertex))
    { }
};

Vertex internal::BasicNonterminalVertex::negate() const {
  return NonterminalVertex(::std::make_shared<BasicNonterminalVertex>(
      ::std::shared_ptr<NonterminalVertexNodeData>(node_data_ptr_),
      !negate_node_
  ));
}

}//namespace slp
}//namespace crag

namespace std {
  //! Definition of the hash for std::pair
  template<typename TFirst, typename TSecond>
  struct hash<pair<TFirst, TSecond>> {
  private:
    const hash<TFirst> first_hash_;
    const hash<TSecond> second_hash_;
  public:
    hash()
      : first_hash_()
      , second_hash_()
    { }
    size_t operator()(const pair<TFirst, TSecond>& obj) const {
      size_t first_hash_value = first_hash_(obj.first);
      //Taken from boost/functional/hash
      return second_hash_(obj.second) + 0x9e3779b9 + (first_hash_value << 6) + (first_hash_value >> 2);
    }
  };

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
