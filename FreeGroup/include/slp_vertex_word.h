/*
 * slp_produced_word.h
 *
 *  Created on: Feb 21, 2013
 *      Author: dpantele
 */

#pragma once
#ifndef CRAG_FREEEGROUP_SLP_VERTEX_WORD_H_
#define CRAG_FREEEGROUP_SLP_VERTEX_WORD_H_

#include "slp_vertex.h"
#include "slp_inspector.h"
#include "slp_pattern_matching.h"

namespace crag {
namespace slp {
namespace inspector {

template <typename AcceptFunctor>
class TerminalVertexPath : public InspectorPath<AcceptFunctor> {
  public:
    TerminalVertexPath() {}
    TerminalVertexPath(AcceptFunctor&& accept_functor)
      : InspectorPath<AcceptFunctor>(std::move(accept_functor))
    { }
    using InspectorPath<AcceptFunctor>::schedule;
    using InspectorPath<AcceptFunctor>::pop_scheduled;
    using InspectorPath<AcceptFunctor>::is_task_accepted;

    InspectorTask process(const InspectorTask& current_task) {
      if (!current_task || current_task.command == InspectorTask::Command::VISIT) {
        return InspectorTask();
      }

      if (current_task.vertex.height() == 1) {
        return InspectorTask::for_current(current_task, InspectorTask::Command::VISIT);
      }

      schedule(InspectorTask::for_right_child(current_task, InspectorTask::Command::GO_LEFT));
      return InspectorTask::for_left_child(current_task, InspectorTask::Command::GO_LEFT);
    }

    InspectorTask initial_task(Vertex&& root) {
      Vertex current = std::move(root);

      while (current.height() > 1) {
        if (current.left_child()) {
          schedule(InspectorTask(current.right_child(), InspectorTask::Command::GO_LEFT, current.left_child().length()));
          current = current.left_child();
        } else {
          current = current.right_child();
        }
      }

      return InspectorTask(std::move(current), InspectorTask::Command::VISIT, LongInteger());
    }

    InspectorTask go_to_position(InspectorTask&& current_task, const LongInteger& new_position) {
      current_task = pop_scheduled();

      while (current_task.vertex.height() > 1 ||
             new_position != current_task.left_siblings_length) {
        if (current_task.left_siblings_length + current_task.vertex.length() < new_position) {
          current_task = pop_scheduled();
        } else {
          current_task = process(current_task);

          if (!is_task_accepted(current_task)) {
            current_task = pop_scheduled();
          }
        }
      }

      return InspectorTask::for_current(current_task, InspectorTask::Command::VISIT);
    }
};

template <typename AcceptFunctor = std::nullptr_t>
class TerminalVertexInspector : public Inspector<TerminalVertexPath, AcceptFunctor> {
  public:
    TerminalVertexInspector(Vertex root)
      : Inspector<TerminalVertexPath, AcceptFunctor>(std::move(root))
    { }

    TerminalVertexInspector() {}
    TerminalVertexInspector(Vertex root, AcceptFunctor accept_functor)
      : Inspector<TerminalVertexPath, AcceptFunctor>(std::move(root), std::move(accept_functor))
    { }

    void go_to_position(const LongInteger& offset) {
      if (offset <= 0) {
        return;
      }

      this->current_task_ = this->inspector_path_.go_to_position(::std::move(this->current_task_), this->current_task_.left_siblings_length + offset);
    }
};

} //namespace inspector

//! Iterator of #slp::VertexWord.
/**
 * This is a forward iterator for slp::VertexWord. Note that even difference_type is standard ptrdiff_t,
 * it is possible that the actual difference between two iterator may be more than 2^sizeof(ptrdiff_t). However,
 * in all standard algorithms this distance can not be achieved, because it takes too long to move
 * this iterator so far.
 */
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
      , current_symbol_()
    { }

    explicit VertexWordIterator(const Vertex& root)
      : inspector_(root)
      , root_(&root)
      , current_symbol_(TerminalVertex(inspector_.vertex()).terminal_symbol())
    { }

    //use default copy/move constructors/assignments

    VertexWordIterator& operator++() {   //!< Preincrement
      ++inspector_;
      current_symbol_ = TerminalVertex(inspector_.vertex()).terminal_symbol();
      return *this;
    }

    VertexWordIterator operator++(int) { //!< Postincrement
      VertexWordIterator copy(*this);
      ++(*this);
      return copy;
    }

    reference operator*() const { //!< "Dereference" current symbol
      return current_symbol_;
    }

    pointer operator->() const {
      return &current_symbol_;
    }

    VertexWordIterator& operator+=(const LongInteger& shift) {
      inspector_.go_to_position(shift);
      current_symbol_ = TerminalVertex(inspector_.vertex()).terminal_symbol();
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
    inspector::TerminalVertexInspector<> inspector_; //!< Subtree inspector
    const Vertex* root_;                           //!< Root. Used only for comparison
    TerminalSymbol current_symbol_;                //!< We have to store it due to type conversion
};

//! Word produced by some #slp::Vertex
class VertexWord {
  public:
    typedef TerminalSymbol value_type;
    typedef VertexWordIterator const_iterator; //no iterator, only const
    typedef const TerminalSymbol& const_reference;
    //we do not define size_type, because size() should return LongInteger, which is not integral

    VertexWord() //!< Just empty word
      : root_(Vertex())
    { }

    template<typename VertexReference>
    explicit VertexWord(VertexReference&& root) //!< Word produced by some root
      : root_(std::forward<VertexReference>(root))
    { }

   //TODO: define swap

    bool operator==(const VertexWord& other) const {
      MatchingTable table;
      return is_equal_to(other, &table);
    }

    bool operator!=(const VertexWord& other) {
       return !(*this == other);
    }

    bool is_equal_to(const VertexWord& other, MatchingTable* matching_table) const {
      if (this->root_.length() != other.root_.length()) {
        return false;
      }

      return static_cast<bool>(matching_table->matches(root_, other.root_));
    }

    inline TerminalSymbol operator[](LongInteger index) const; //!< Get one letter from the word

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

inline ::std::ostream& operator<<(::std::ostream& out, const VertexWord& word) {
  for(auto symbol : word) {
    out << symbol;
  }

  return out;
}

inline TerminalSymbol VertexWord::operator[](LongInteger index) const {
  Vertex current_vertex = root_;

  while (current_vertex.height() > 1) {
    if (index < current_vertex.left_child().length()) { //Go left
      current_vertex = current_vertex.left_child();
    } else {
      index -= current_vertex.left_child().length();
      current_vertex = current_vertex.right_child();
    }
  }

  return TerminalVertex(current_vertex).terminal_symbol();
}

} // namespace slp
} // namespace crag

#endif /* SLP_VERTEX_WORD_H_ */
