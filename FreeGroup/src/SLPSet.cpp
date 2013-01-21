#include "SLPSet.h"

const SLPVertex SLPVertex::Null; //define static in cpp file

SLPVertex SLPVertex::terminal_vertex(const TerminalSymbol& terminal_symbol) {
  SLPVertex vertex(std::make_shared<BasicVertex>(), false);

  vertex.ptr_->terminal_symbol = terminal_symbol;
  vertex.ptr_->length = 1;
  vertex.ptr_->height = 1;

  return vertex;
}

/* TODO: implement equal_to
SLPSet::equal_to(const SLPSet& other) const {
  if (other.roots.size() != this->roots.size()) {
    return false;
  }

  if (other.terminals_count != this->terminals_count) {
    return false;
  }

  for (auto this_root = this->roots.begin(), auto other_root = other.roots.end();
       this_root != this->roots.end(); //Don't have to check other, because other.roots.size() == this->roots.size()
       ++this_root, ++other_root) {
    if (this->vertices[*this_root].length != other.vertices[*other_root].length) {
      return false;
    }
  }

  SLPMatchingTable PT(*this, other);
  for (auto this_root = this->roots.begin(), auto other_root = other.roots.end();
       this_root != this->roots.end(); //Don't have to check other, because other.roots.size() == this->roots.size()
       ++this_root, ++other_root) {
     
    SLPMatchingTable::MatchResultSequence match = PT.get_matches(SignedVertex(*this_root, false), SignedVertex(*other_root, false));
    if (match.start != 0 || match.count != 1) {
      return false;
    }

  }

  return true;

}
*/
