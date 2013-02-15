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

//! Basic interface for vertex. Also represents empty vertex, use Vertex::Null as empty vertex
class Vertex {
  public:
    typedef ::std::shared_ptr<Vertex> Ptr;

    virtual ~Vertex() {} //We are required to have virtual destructor

    const static Ptr Null; //!< Use this vertex to represent invalid (or empty) vertex

    bool operator==(const Vertex& other) const {
      return typeid(*this) == typeid(other) && call_is_same_vertex(other);
    }

    bool operator!=(const Vertex& other) const {
      return !(*this == other);
    }

    //! Return
    virtual Ptr negate() const {
      return Null;
    }

    virtual Ptr left_child() const {
      return Null;
    }

    virtual bool has_left_child() const {
      return false;
    }

    virtual Ptr right_child() const {
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

    virtual Ptr clone() const {
      return Null;
    }

    static Ptr create() {
      if (Null) {
        return Null;
      }
      return ::std::shared_ptr<Vertex>(new Vertex());
    }

  protected:
    Vertex() {} //!< Default constructor generating empty vertex. Use #Null instead of it

    virtual bool call_is_same_vertex(const Vertex& other) const {
      return true;
    }

    bool is_same_vertex(const Vertex& other) const {
      return true;
    }
};

const Vertex::Ptr Vertex::Null = Vertex::create();

//! Terminal vertex in a SLP. Produces word of length 1.
template <typename TerminalSymbol>
class TerminalVertexTemplate : public Vertex {
  public:
    TerminalVertexTemplate() = delete;

    static Ptr create(const TerminalSymbol& terminal_symbol) {
      return ::std::shared_ptr<TerminalVertexTemplate>(new TerminalVertexTemplate(terminal_symbol));
    }

    virtual Ptr negate() const {
      return ::std::shared_ptr<TerminalVertexTemplate>(new TerminalVertexTemplate(-terminal_symbol_));
    }

    virtual LongInteger length() const {
      return 1;
    }

    virtual unsigned int height() const {
      return 1;
    }

    virtual size_t vertex_hash() const {
      return hash_(terminal_symbol_);
    }

    const TerminalSymbol& terminal_symbol() const {
      return terminal_symbol();
    }

    operator TerminalSymbol() const {
      return terminal_symbol();
    }

    virtual Ptr clone() const {
      return ::std::make_shared<TerminalVertexTemplate>(*this);
    }

  protected:

    virtual bool call_is_same_vertex(const Vertex& other) const {
      return is_same_vertex(dynamic_cast<const TerminalVertexTemplate&>(other));
    }

    bool is_same_vertex(const TerminalVertexTemplate& other) const {
      return this->terminal_symbol_ == other.terminal_symbol_;
    }
  private:
    ::std::hash<TerminalSymbol> hash_;
    TerminalSymbol terminal_symbol_;

    explicit TerminalVertexTemplate(const TerminalSymbol& terminal_symbol)
      : terminal_symbol_(terminal_symbol)
    { }

};

template <typename TerminalSymbol>
::std::ostream& operator << (::std::ostream& stream, const TerminalVertexTemplate<TerminalSymbol>& vertex) {
  return stream << vertex.terminal_symbol();
}

//! Non-terminal vertex in a SLP, represent rule A->BC
class NonterminalVertex : public Vertex {
  private:
    struct BasicVertex {
      Ptr left_child;            //!< Left part of the rule. Use SignedVertex::Null if child is absent.
      Ptr right_child;           //!< Right part of the rule.
      LongInteger length;              //!< Length of word produced by the vertex
      unsigned int height;             //!< Height of subtree

      BasicVertex(Ptr left_child, Ptr right_child, LongInteger length, unsigned int height)
        : left_child(left_child)
        , right_child(right_child)
        , length(length)
        , height(height)
      { }
    };
  public:
    NonterminalVertex() = delete;

    static Ptr create(const Ptr& left_child, const Ptr& right_child) {
      return ::std::shared_ptr<NonterminalVertex>(new NonterminalVertex(
          ::std::make_shared<BasicVertex>(
              left_child->clone(),
              right_child->clone(),
              left_child->length() + right_child->length(),
              ::std::max(left_child->height(), right_child->height()) + 1
          ),
          false));
    }

    virtual Ptr negate() const {
      return ::std::shared_ptr<NonterminalVertex>(new NonterminalVertex(basic_vertex_, !is_negative_));
    }

    virtual Ptr left_child() const {
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

    virtual Ptr right_child() const {
      if (is_negative_) {
        return basic_vertex_->left_child->negate();
      } else {
        return basic_vertex_->right_child;
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

    virtual size_t vertex_hash() const {
      return is_negative_ ? (~ptr_hash_(basic_vertex_)) : ptr_hash_(basic_vertex_);
    }

    virtual Ptr clone() const {
      return ::std::make_shared<NonterminalVertex>(*this);
    }

  protected:

    virtual bool call_is_same_vertex(const Vertex& other) const {
      return is_same_vertex(dynamic_cast<const NonterminalVertex&>(other));
    }

    bool is_same_vertex(const NonterminalVertex& other) const {
      return this->basic_vertex_ == other.basic_vertex_ && this->is_negative_ == other.is_negative_;
    }

  private:
    ::std::hash<std::shared_ptr<BasicVertex>> ptr_hash_;
    ::std::shared_ptr<BasicVertex> basic_vertex_;
    bool is_negative_;

    NonterminalVertex(::std::shared_ptr<BasicVertex> basic_vertex, bool is_negative)
      : basic_vertex_(basic_vertex)
      , is_negative_(is_negative)
    {  }

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


#endif /* CRAG_FREEGROUP_SLP_VERTEX_H_ */
