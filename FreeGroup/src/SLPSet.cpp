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

const SLPMatchingTable::MatchResultSequence SLPMatchingTable::NO_MATCHES;

bool terminals_equal(const SLPVertex& lhs, const SLPVertex& rhs) {
  return lhs.terminal_symbol() == rhs.terminal_symbol() && lhs.is_negative() == rhs.is_negative();
}

SLPMatchingTable::MatchResultSequence SLPMatchingTable::matches(const SLPVertex& pattern,
                                              const SLPVertex& text) {
  auto match_result_iterator = match_table_.find(std::make_pair(pattern, text));

  if (match_result_iterator != match_table_.end()) { //if already calculated
    return match_result_iterator->second;
  }

  if (pattern.length() > text.length()) {
    return NO_MATCHES;
  }

  MatchResultSequence match_result;

  if (pattern.length() == 1) {//Trivial case
    SLPVertex pattern_vertex = SLPProducedWord(pattern)[0];
    if (text.length() == 1) {
      if (terminals_equal(pattern_vertex, SLPProducedWord(text)[0])) {
        match_result.count = 1;
        match_result.start = 0;
        match_result.step = 1;
      } else {
        match_result = NO_MATCHES;
      }
    } else {
      SLPProducedWord left_text_part = SLPProducedWord(text.left_child());
      SLPProducedWord right_text_part = SLPProducedWord(text.right_child());

      bool left_match = left_text_part.size() > 0 &&
          terminals_equal(pattern_vertex, left_text_part[left_text_part.size() - 1]);
      bool right_match = right_text_part.size() > 0 &&
          terminals_equal(pattern_vertex, right_text_part[0]);

      if (left_match) {
        match_result.count = right_match ? 2 : 1;
        match_result.start = left_text_part.size() - 1;
        match_result.step = 1;
      } else if (right_match) {
        match_result.count = 1;
        match_result.start = left_text_part.size();
        match_result.step = 1;
      } else {
        match_result = NO_MATCHES;
      }
    }
  } else  {//we have pattern.length > 1 => text.length > 1
    if (pattern.left_child().length() >= pattern.right_child().length()) {//Right child is smaller
      auto local_matches = local_search(pattern.left_child(), text,
                                        text.left_child().length() - pattern.length(),
                                        text.left_child().length() + pattern.left_child().length());
      //TODO: finish this mess
    }
  }

  match_table_.insert(std::make_pair(std::make_pair(pattern, text), match_result));

  return match_result;
}

class SLPMatchingInspector {
  public:
    SLPMatchingInspector() = delete;
    SLPMatchingInspector(SLPVertex pattern, SLPVertex text, bool reverse, LongInteger left_border, LongInteger right_border)
      : parent_stack_()
      , pattern_(pattern)
      , text_position_({text, 0})
      , reverse_(reverse)
      , left_border_(left_border)
      , right_border_(right_border)
    {
      go_further();
    }

    bool inspection_ended() const {
      return text_position_.vertex == SLPVertex::Null && parent_stack_.empty();
    }

    void go_next() {
      if (reverse_) {
        text_position_.vertex = text_position_.vertex.left_child();
      } else {
        text_position_.text_begin += text_position_.vertex.left_child().length();
        text_position_.vertex = text_position_.vertex.right_child();
      }
      go_further();
    }
  private:
    struct TextPosition {
      SLPVertex vertex;
      LongInteger text_begin;
    };
    std::vector<TextPosition> parent_stack_;
    SLPVertex pattern_;
    TextPosition text_position_;
    bool reverse_;
    LongInteger left_border_;
    LongInteger right_border_;

    bool position_is_suitable(const TextPosition& position, bool is_right_child) const {
      if (text_position_.vertex == SLPVertex::Null ||
          right_border_ - text_position_.text_begin < pattern_.length() ||
          text_position_.text_begin + text_position_.vertex.length() - left_border_ < pattern_.length()
          ) {
        return false;
      }

      return true;
    }

    void go_further() {
      //TODO: case root.length() == pattern_.length()
      while (position_is_suitable(text_position_, reverse_)) {
        parent_stack_.push_back(text_position_);

        if (reverse_) {
          text_position_.text_begin += text_position_.vertex.left_child().length();
          text_position_.vertex = text_position_.vertex.right_child();
        } else {
          text_position_.vertex = text_position_.vertex.left_child();
        }
      }
      if (!parent_stack_.empty()) {
        text_position_ = parent_stack_.back();
        parent_stack_.pop_back();
      } else {
        text_position_ = {SLPVertex::Null, 0};
      }
    }

};

} //namespace crag