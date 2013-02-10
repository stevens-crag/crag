/*
 * slp_vertex.h
 *
 *  Created on: Feb 10, 2013
 *      Author: dpantele
 */

#ifndef SLP_VERTEX_H_
#define SLP_VERTEX_H_

#include <gmpxx.h>
typedef mpz_class LongInteger;

#include <memory>

namespace crag {
namespace slp {

//! Basic interface for vertex. Also represents empty vertex, use Vertex::Null as empty vertex
class Vertex {
  public:
    typedef ::std::shared_ptr<Vertex> VertexPtr;
    typedef ::std::shared_ptr<const Vertex> ConstVertexPtr;

    virtual ~Vertex() {} //We are required to have virtual destructor

    const static VertexPtr Null; //!< Use this vertex to represent invalid (or empty) vertex

    bool operator==(const Vertex& other) const {
      return typeid(*this) == typeid(other) && call_is_same_vertex(other);
    }

    bool operator!=(const Vertex& other) const {
      return !(*this == other);
    }

    //! Return
    virtual VertexPtr negate() const {
      return Null;
    }

    virtual VertexPtr left_child() const {
      return Null;
    }

    virtual bool has_left_child() const {
      return false;
    }

    virtual VertexPtr right_child() const {
      return Null;
    }

    virtual bool has_right_child() const {
      return false;
    }

    virtual LongInteger length() const {
      return 0;
    }

    virtual unsigned int height() const {
      return 0;
    }

    virtual size_t vertex_hash() const {
      return 0;
    }

  private:
    Vertex() {} //!< Default constructor generating empty vertex. Use #Null instead of it

    virtual VertexPtr clone() const {
      return Null;
    }

    virtual bool call_is_same_vertex(const Vertex& other) const {
      return true;
    }

    bool is_same_vertex(const Vertex& other) {
      return true;
    }
};

const static Vertex::VertexPtr Vertex::Null = ::std::make_shared<Vertex>();

//! Terminal vertex in a SLP. Produces word of length 1.
template <typename TerminalSymbol>
class TerminalVertex : public Vertex {
  public:
    TerminalVertex() = delete;

    static VertexPtr create(const TerminalSymbol& terminal_symbol) {
      return ::std::make_shared<TerminalVertex>(terminal_symbol);
    }

    virtual Vertex negate() const {
      return ::std::make_shared<TerminalVertex>(-terminal_symbol_);
    }

    virtual LongInteger length() const {
      return 1;
    }

    virtual unsigned int height() const {
      return 1;
    }

    virtual size_t vertex_hash() {
      return hash_(terminal_symbol_);
    }

    const TerminalSymbol& terminal_symbol() const {
      return terminal_symbol();
    }

    operator TerminalSymbol() const {
      return terminal_symbol();
    }

  private:
    ::std::hash<TerminalSymbol> hash_;
    TerminalSymbol terminal_symbol_;

    explicit TerminalVertex(const TerminalSymbol& terminal_symbol)
      : terminal_symbol_(terminal_symbol)
    { }

    virtual VertexPtr clone() const {
      return ::std::make_shared<TerminalVertex>(*this);
    }

    virtual bool call_is_same_vertex(const Vertex& other) const {
      return is_same_vertex(dynamic_cast<const TerminalVertex&>(other));
    }

    bool is_same_vertex(const TerminalVertex& other) {
      return this->terminal_symbol_ == other.terminal_symbol_;
    }
};

template <typename TerminalSymbol>
::std::ostream& operator << (::std::ostream& stream, const TerminalVertex<TerminalSymbol>& vertex) {
  return ::std::ostream << vertex.terminal_symbol();
}

//! Non-terminal vertex in a SLP, represent rule A->BC
class NonterminalVertex : public Vertex {
  private:
    struct BasicVertex {
      VertexPtr left_child;            //!< Left part of the rule. Use SignedVertex::Null if child is absent.
      VertexPtr right_child;           //!< Right part of the rule.
      LongInteger length;              //!< Length of word produced by the vertex
      unsigned int height;             //!< Height of subtree

      BasicVertex(VertexPtr left_child, VertexPtr right_child, LongInteger length, LongInteger height)
        : left_child(left_child)
        , right_child(right_child)
        , length(length)
        , height(height)
      { }
    };
  public:
    NonterminalVertex() = delete;

    static VertexPtr create(const VertexPtr& left_child, const VertexPtr& right_child) {
      return ::std::make_shared<NonterminalVertex>(
          BasicVertex(
              left_child->clone(),
              right_child->clone(),
              left_child->length() + right_child->length(),
              ::std::max(left_child->height(), right_child->height()) + 1
          ),
          false);
    }

    virtual VertexPtr negate() const {
      return ::std::make_shared<NonterminalVertex>(basic_vertex_, !is_negative_);
    }

    virtual VertexPtr left_child() const {
      if (is_negative_) {
        return basic_vertex_->right_child->negate();
      } else {
        return basic_vertex_->left_child;
      }
    }

    virtual bool has_left_child() const {
      if (is_negative_) {
        return basic_vertex_->right_child != Null;
      } else {
        return basic_vertex_->left_child != Null;
      }
    }

    virtual VertexPtr right_child() const {
      if (is_negative_) {
        return basic_vertex_->right_child->negate();
      } else {
        return basic_vertex_->left_child;
      }
    }

    virtual bool has_right_child() const {
      if (is_negative_) {
        return basic_vertex_->left_child != Null;
      } else {
        return basic_vertex_->right_child != Null;
      }
    }

    virtual LongInteger length() const {
      return basic_vertex_->length;
    }

    virtual unsigned int height() const {
      return basic_vertex_->height;
    }

    virtual size_t vertex_hash() {
      return is_negative_ ? (~ptr_hash_(basic_vertex_)) : ptr_hash_(basic_vertex_);
    }

  private:
    ::std::hash<std::shared_ptr<BasicVertex>> ptr_hash_;
    ::std::shared_ptr<BasicVertex> basic_vertex_;
    bool is_negative_;

    NonterminalVertex(::std::shared_ptr<BasicVertex> basic_vertex, bool is_negative)
      : basic_vertex_(basic_vertex)
      , is_negative_(is_negative)
    {  }

    virtual VertexPtr clone() const {
      return ::std::make_shared<TerminalVertex>(*this);
    }

    virtual bool call_is_same_vertex(const Vertex& other) const {
      return is_same_vertex(dynamic_cast<const NonterminalVertex&>(other));
    }

    bool is_same_vertex(const NonterminalVertex& other) {
      return this->basic_vertex_ == other.basic_vertex_ && this->is_negative_ && other.is_negative_;
    }
};
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


#endif /* SLP_VERTEX_H_ */
