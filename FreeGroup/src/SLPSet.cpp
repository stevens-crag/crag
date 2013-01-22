#include "SLPSet.h"

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
