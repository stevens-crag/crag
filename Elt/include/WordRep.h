// Copyright (C) 2005 Alexander Ushakov
// Contents: Definition of class WordRep
//
// Principal Authors: Alexander Ushakov
//

#ifndef CRAG_WORDREP_H
#define CRAG_WORDREP_H

#include <algorithm>
#include <list>
#include <ostream>
#include <stdexcept>
#include <vector>

//! Represents a group word.
//! Usually a word is kept reduced, except for the cases when insert/replace functions are used.
class WordRep {
private:
  using storage_t = std::vector<int>;

public:
  using iterator = storage_t::iterator;
  using const_iterator = storage_t::const_iterator;
  using const_reverse_iterator = storage_t::const_reverse_iterator;

  WordRep() {}

  explicit WordRep(int g);

  explicit WordRep(const std::list<int>& gens);

  explicit WordRep(std::vector<int> gens);

  explicit WordRep(std::initializer_list<int> gens);

  template<class InputIterator>
  WordRep(InputIterator begin, InputIterator end);

  std::string toString() const;

  iterator begin();

  const_iterator begin() const;

  iterator end();

  const_iterator end() const;

  const_iterator cbegin() const;

  const_iterator cend() const;

  const_reverse_iterator rbegin() const;

  const_reverse_iterator rend() const;

  int front() const {
    return elements_.front();
  }

  int& front() {
    return elements_.front();
  }

  int back() const {
    return elements_.back();
  }

  int& back() {
    return elements_.back();
  }

  // relational operators
  bool operator==(const WordRep& other) const;

  bool operator!=(const WordRep& other) const;

  bool operator<(const WordRep& other) const;

  bool operator>(const WordRep& other) const;

  //! Conjugate a word by another word
  WordRep& operator^=(const WordRep& conjugator);

  inline WordRep operator^(const WordRep& conjugator) const {
    WordRep result(*this);
    result ^= conjugator;
    return result;
  }

  //! Raise a word into a power
  WordRep& operator^=(int power);

  inline WordRep operator^(int power) const {
    WordRep result(*this);
    result ^= power;
    return result;
  }

  WordRep& operator*=(const WordRep& other);

  inline WordRep operator*(const WordRep& w) const {
    WordRep result(*this);
    result *= w;
    return result;
  }

  bool contains(int gen) const;

  inline size_t length() const {
    return elements_.size();
  }

  inline size_t size() const {
    return elements_.size();
  }

  inline bool empty() const {
    return elements_.empty();
  }

  //! Counts exponent sum of x_{|gen|}
  int exponentSum(int gen) const;

  //! Counts occurrences of x_{|gen|}^{+/- 1} in this word.
  size_t occurrences(int gen) const;

  //! Extracts root from this word.
  std::pair<WordRep, int> root() const;

  WordRep inverse() const;

  //! Make a word trivial
  inline void clear() {
    elements_.clear();
  }

  //! Freely reduces this word.
  void freelyReduce();

  //! Freely reduces the fragment [begin, end) of this word.
  void freelyReduce(iterator begin, iterator end);

  //! Cyclically reduces this word w to w' and returns g such that g^{-1} w' g = w
  WordRep cyclicallyReduce();

  //! the leftmost symbol goes to the rightmost position
  void cyclicLeftShift();

  //! the rightmost symbol goes to the leftmost position
  void cyclicRightShift();

  //! First makes cyclic permutations and only then makes free reduction.
  //! So it is not equal to a composition of several left/right cyclic shifts.
  //! If n > 0 => left-shift permute
  //! If n < 0 => rigth-shift permute
  void cyclicallyPermute(int n);

  //! Replaces this word with its segment [from, to).
  void segment(size_t from, size_t to);

  //! Returns subword [from, to).
  WordRep subword(size_t from, size_t to) const;

  //! Replaces this word with its prefix [0, to).
  void initialSegment(size_t to);

  //! Replaces this word with its suffix [from, size - 1).
  void terminalSegment(size_t from);

  //! Inserts a range of nonzero elements at the given position.
  //! Reduction is not performed!
  template<class InputIterator>
  void insert(size_t position, InputIterator begin, InputIterator end);

  //! Inserts a range of nonzero elements at the given position.
  //! Reduction is not performed!
  template<class InputIterator>
  void insert(iterator position, InputIterator begin, InputIterator end);

  //! Inserts element g \neq 0 at the given position.
  //! Inserts at the end if pos >= size().
  //! Reduction is not performed!
  void insert(size_t position, int g);

  //! Inserts element g \neq 0 at the given position.
  //! Reduction is not performed!
  void insert(iterator position, int g);

  //! Replaces element at the given position with g \neq 0.
  //! Throws if pos >= size().
  //! Reduction is not performed!
  void replace(size_t position, int g);

  //! Replaces element at the given position with g \neq 0.
  //! Reduction is not performed!
  void replace(iterator position, int g);

  //! Replaces a subword of a word starting at a position it by a word [b, e).
  //! The length of the word does not increase if [b, e) is longer than the terminal segment of the word [it, end()).
  //! In that case terminal symbols of [b, e) are ignored.
  //! Reduction is not performed!
  template<class InputIterator>
  void replace(iterator position, InputIterator b, InputIterator e);

  template<class InputIterator>
  void replace(size_t position, InputIterator b, InputIterator end);

  void push_back(int g) {
    reduced_push_back_(g);
  }

  void pop_back() {
    elements_.pop_back();
  }

  void push_front(int g) {
    reduced_push_front_(g);
  }

  void pop_front() {
    elements_.erase(elements_.begin());
  }

  std::list<int> toList() const {
    return std::list<int>(elements_.begin(), elements_.end());
  }

  const std::vector<int>& toVector() const {
    return elements_;
  }

private:
  inline void validate_(int g) const {
    if (g == 0) {
      throw std::invalid_argument("Zero indices are not allowed.");
    }
  }

  template<typename Iterator>
  void validate_(Iterator begin, Iterator end) const {
    for (auto it = begin; it != end; ++it) {
      validate_(*it);
    }
  }

  void validate_() const {
    validate_(elements_.begin(), elements_.end());
  }

  //! Multiply on the right and reduce
  void reduced_push_back_(int g);

  //! Multiply on the left and reduce
  void reduced_push_front_(int g);

  // list of generators, negative integers represent inverses of positive integers
  storage_t elements_;
};

template<class InputIterator>
WordRep::WordRep(InputIterator begin, InputIterator end){
  validate_(begin, end);

  for (auto it = begin; it != end; ++it) {
    reduced_push_back_(*it);
  }
}

template<class InputIterator>
void WordRep::insert(size_t position, InputIterator begin, InputIterator end) {
  position = std::min(position, size());

  auto it = elements_.begin();
  std::advance(it, position);

  insert(it, begin, end);
}

template<class InputIterator>
void WordRep::insert(iterator position, InputIterator begin, InputIterator end) {
  validate_(begin, end);

  elements_.insert(position, begin, end);
}

template<class InputIterator>
void WordRep::replace(iterator position, InputIterator b, InputIterator e) {
  for (; position != end() && b != e; ++position, ++b) {
    *position = *b;
  }
}

template<class InputIterator>
void WordRep::replace(size_t position, InputIterator b, InputIterator e) {
  if (position >= size()) {
    throw std::invalid_argument("Bad position.");
  }

  auto it = elements_.begin();
  std::advance(it, position);

  replace(it, b, e);
}

std::ostream& operator<<(std::ostream& out, const WordRep& w);

#endif // CRAG_WORDREP_H
