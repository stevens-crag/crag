#include "SLPSet.h"

namespace crag {

const SLPVertex SLPVertex::Null; //define static in cpp file

inline BasicVertex::BasicVertex()
  : left_child(SLPVertex::Null),
    right_child(SLPVertex::Null),
    terminal_symbol(INVALID_TERMINAL),
    length(0),
    height(0)
{ }

SLPVertex SLPVertex::terminal_vertex(const TerminalSymbol& terminal_symbol) {
  SLPVertex vertex(std::make_shared<BasicVertex>(), false);

  vertex.ptr_->terminal_symbol = terminal_symbol;
  vertex.ptr_->length = 1;
  vertex.ptr_->height = 1;

  return vertex;
}

SLPVertex SLPVertex::concatenate(const SLPVertex& left_child, const SLPVertex& right_child) {
  SLPVertex vertex(std::make_shared<BasicVertex>(), false);

  vertex.ptr_->left_child = left_child;
  vertex.ptr_->right_child = right_child;
  vertex.ptr_->length = left_child.length() + right_child.length();
  vertex.ptr_->height = std::max(left_child.height(), right_child.height()) + 1;

  return vertex;
}

void SLPPostorderInspector::go_to_next_vertex() {
  //TODO: maybe raise an exception if current_path is empty already?
  const SLPVertex current = current_path_.back();
  current_path_.pop_back();

  if (inspection_ended()) {
    return;
  }

  const SLPVertex& parent = current_path_.back();

  if ( parent.right_child() == current || !parent.has_right_child()) {
    //We have already visited all right children of the parent vertex, stop on parent
  } else {
    //We have visited all left children of parent, going to the right
    current_path_.push_back(parent.right_child());
    goto_leftmost_terminal();
  }
}

void SLPPostorderInspector::goto_leftmost_terminal() {
  while (current_path_.back() != SLPVertex::Null) {
    const SLPVertex & current = current_path_.back();
    if (current.has_left_child()) {
      current_path_.push_back(current.left_child());
    } else {
      //Don't have any left children, then go right
      current_path_.push_back(current.right_child());
    }
  }

  current_path_.pop_back();
}

SLPProducedWordIterator& SLPProducedWordIterator::operator ++() {
  ++length_;
  inspector_.go_to_next_vertex();

  while (!inspector_.inspection_ended() && !inspector_.current_vertex().is_terminal()) {
    inspector_.go_to_next_vertex();
  }

  return *this;
}

SLPProducedWord::value_type SLPProducedWord::operator [](LongInteger index) const {
  SLPVertex current_vertex = root_;

  while (!current_vertex.is_terminal()) {
    if (index < current_vertex.left_child().length()) { //Go left
      current_vertex = current_vertex.left_child();
    } else {
      index -= current_vertex.left_child().length();
      current_vertex = current_vertex.right_child();
    }
  }

  return current_vertex;
}

} //namespace crag
