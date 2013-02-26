/*
 * slp_produced_word.h
 *
 *  Created on: Feb 21, 2013
 *      Author: dpantele
 */

#ifndef CRAG_FREEEGROUP_SLP_VERTEX_WORD_H_
#define CRAG_FREEEGROUP_SLP_VERTEX_WORD_H_

#include "slp_vertex.h"
#include "slp_inspector.h"

namespace crag {
namespace slp {

namespace internal {

template <typename OrderPolicy>
class TerminalVertexInspector : public Inspector<OrderPolicy> {
    typedef Inspector<OrderPolicy> Super;
  public:
    template<typename ... ConstructorArgumentTypes>
    TerminalVertexInspector(ConstructorArgumentTypes&&... arguments)
      : Super(std::forward<ConstructorArgumentTypes>(arguments)...)
    {
      if (this->current_task_.vertex.height() > 1) {
        next();
      }
    }

    TerminalVertexInspector& next() {
      do {
        Super::next();
      } while (this->current_task_.vertex.height() > 1); //Null vertex has height 0, so if current_task_ == DO_NOTHING, this loop breaks
      return *this;
    }

    void go_to_position(const LongInteger& position) {
      if (position < this->current_task_.left_siblings_length) {
        return; //We have to go backwards
      }

      while (this->current_task_.vertex.height() > 1 ||
             this->current_task_.command != Super::InspectorTask::Command::VISIT ||
             position != this->current_task_.left_siblings_length) {
        if (this->current_task_.left_siblings_length + this->current_task_.vertex.length() < position) {
          this->fallback_to_scheduled(); //We have to skip current subtree
        } else {
          this->process_current_task();
        }
      }
      //While we are not visiting terminal or null
    }
};

} //namespace internal

//! Iterator of #slp::VertexWord.
/**
 * This is a forward iterator for slp::VertexWord. Note that even difference_type is standard ptrdiff_t,
 * it is possible that the actual difference between two iterator may be more than 2^sizeof(ptrdiff_t). However,
 * in all standard algorithms this distance can not be achieved, because it takes too long to move
 * this iterator so far.
 */
template <typename TerminalSymbol>
class VertexWordIterator : public std::iterator <
        std::forward_iterator_tag,      //iterator_category
        const TerminalSymbol            //value_type
  > {
  public:
    typedef std::forward_iterator_tag  iterator_category;
    typedef const TerminalSymbol       value_type;
    typedef ::std::ptrdiff_t           difference_type;
    typedef const TerminalSymbol*      pointer;
    typedef const TerminalSymbol&      reference;

    VertexWordIterator()
      : inspector_()
      , root_(nullptr)
    { }

    explicit VertexWordIterator(const Vertex& root)
      : inspector_(root)
      , root_(&root)
    { }

    //use default copy/move constructors/assignments

    VertexWordIterator& operator++() {   //!< Preincrement
      ++inspector_;
      return *this;
    }

    VertexWordIterator operator++(int) { //!< Postincrement
      VertexWordIterator copy(*this);
      ++(*this);
      return copy;
    }

    reference operator*() const { //!< "Dereference" current symbol
      return TerminalVertexTemplate<TerminalSymbol>(inspector_.vertex()).terminal_symbol();
    }

    pointer operator->() const {
      return &TerminalVertexTemplate<TerminalSymbol>(inspector_.vertex()).terminal_symbol();
    }

    VertexWordIterator& operator+=(const LongInteger& shift) {
      inspector_.go_to_position(inspector_.vertex_left_siblings_length() + shift);
      return *this;
    }

    VertexWordIterator operator+(const LongInteger& shift) {
      return VertexWordIterator(*this) += shift;
    }

    //!< Compare to another iterator.
    /**
     * Compares position of inspectors && root. If this->length >= root.length() and other->length >= root.length(), then also true
     */
    bool operator==(const VertexWordIterator& other) const {
      return ((!root_ && other.inspector_.stopped()) || (!other.root_ && inspector_.stopped())) ||
          ( root_ && other.root_ && *root_ == *other.root_ && inspector_ == other.inspector_);
    }

    bool operator!=(const VertexWordIterator& other) const {
      return !(*this == other);
    }

  private:
    internal::TerminalVertexInspector<inspector::Preorder> inspector_; //!< Subtree inspector
    const Vertex* root_;                                               //!< Root. Used only for comparison
};

//! Word produced by some #slp::Vertex
template <typename TerminalSymbol>
class VertexWord {
  public:
    typedef TerminalSymbol value_type;
    typedef VertexWordIterator<TerminalSymbol> const_iterator; //no iterator, only const
    typedef const TerminalSymbol& const_reference;
    //we do not define size_type, because size() should return LongInteger, which is not integral

    VertexWord() //!< Just empty word
      : root_(Vertex::Null)
    { }

    template<typename VertexReference>
    explicit VertexWord(VertexReference&& root) //!< Word produced by some root
      : root_(std::forward<VertexReference>(root))
    { }

   //TODO: define swap, operator ==

    const_reference operator[](LongInteger index) const; //!< Get one letter from the word

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
    Vertex root_; //!< The root vertex producing this word
};

template <typename TerminalSymbol>
typename VertexWord<TerminalSymbol>::const_reference VertexWord<TerminalSymbol>::operator[](LongInteger index) const {
  Vertex current_vertex = root_;

  while (current_vertex.height() > 1) {
    if (index < current_vertex.left_child().length()) { //Go left
      current_vertex = current_vertex.left_child();
    } else {
      index -= current_vertex.left_child().length();
      current_vertex = current_vertex.right_child();
    }
  }

  return TerminalVertexTemplate<TerminalSymbol>(current_vertex).terminal_symbol();
}

} // namespace slp
} // namespace crag

#endif /* SLP_VERTEX_WORD_H_ */
