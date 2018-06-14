// Copyright (C) 2005 Alexander Ushakov
// Contents: Definition of class Word
//
// Principal Authors: Alexander Ushakov
//

#ifndef CRAG_WORD_H
#define CRAG_WORD_H

#include <string>
#include <set>
#include <map>
#include <memory>
#include <tuple>

#include "WordRep.h"
#include "Alphabet.h"

//! Class Word - defines a representation of a Word over a group alphabet.
/*!
  A reduced word is the one which does not involve subwords of the type \f$x x^{-1}\f$ or \f$x^{-1} x\f$.
  We represent a reduced word over a group alphabet \f$X\f$ as a list of non-trivial integers std::list<int>.
  Each generator \f$x \in X\f$ is represented by the unique nonzero number \f$n_x\f$.
  For each \f$x \in X\f$ it is assumed that \f$n_{x^{-1}} = -n_x\f$.
*/
class Word {
private:
  explicit Word(WordRep w)
    : impl_ptr_(std::make_shared<WordRep>(std::move(w))) {}

public:
  using iterator = WordRep::iterator;
  using const_iterator = WordRep::const_iterator;
  using const_reverse_iterator = WordRep::const_reverse_iterator;

  //! Default constructor (creates the empty word \f$w = \varepsilon\f$).
  Word()
    : impl_ptr_(std::make_shared<WordRep>()) {}

  //! Constructs a word by its presentation (if the word defined in gens is not reduced then it reduces it).
  explicit Word(std::vector<int> gens)
    : impl_ptr_(std::make_shared<WordRep>(std::move(gens))) {}

  //! Constructs a word by its presentation (if the word defined in gens is not reduced then it reduces it).
  explicit Word(const std::list<int>& gens)
    : impl_ptr_(std::make_shared<WordRep>(gens)) {}

  explicit Word(std::initializer_list<int> gens)
    : impl_ptr_(std::make_shared<WordRep>(gens)) {}

  //! Constructs a one letter word.
  explicit Word(int g)
    : impl_ptr_(std::make_shared<WordRep>(g)) {}

  template<class InputIterator>
  Word(InputIterator begin, InputIterator end)
    : impl_ptr_(std::make_shared<WordRep>(begin, end)) {}

  //! Comparison operator
  //! The result of comparison of u = u_1...u_k and v = v_1...v_m is defined by
  //! u < v iff k < m or (k = m and u_i < v_i where i is the smallest index such that u_i \neq v_i)
  bool operator<(const Word& other) const {
    return *impl_ptr_ < *other.impl_ptr_;
  }

  //! Comparison operator
  bool operator>(const Word& other) const {
    return other < *this;
  }

  //! Comparison operator
  bool operator==(const Word& other) const {
    return *impl_ptr_ == *other.impl_ptr_;
  }

  //! Comparison operator
  bool operator!=(const Word& other) const {
    return !(*this == other);
  }

  //! Multiply the word on the right by another word (the result is reduced)
  Word& operator*=(const Word& other) {
    clone_();
    *impl_ptr_ *= *other.impl_ptr_;

    return *this;
  }

  //! Multiply two words (the result is reduced)
  Word operator*(const Word& other) const {
    Word result(*this);
    result *= other;
    return result;
  }

  //! Conjugate a word by another word (the result is reduced)
  Word& operator^=(const Word& conjugator) {
    clone_();
    *impl_ptr_ ^= *conjugator.impl_ptr_;

    return *this;
  }

  Word operator^(const Word& conjugator) const {
    Word result(*this);
    result ^= conjugator;
    return result;
  }

  //! Conjugate a word by another word (the result is reduced)
  Word& operator^=(int power) {
    clone_();
    *impl_ptr_ ^= power;

    return *this;
  }

  Word operator^(int power) const {
    Word result(*this);
    result ^= power;
    return result;
  }

  //! Invert a word (works the same as inverse).
  Word operator-() const {
    WordRep inverse = impl_ptr_->inverse();
    return Word(std::move(inverse));
  }

  const_iterator begin() const {
    return impl_ptr_->begin();
  }

  const_iterator end() const {
    return impl_ptr_->end();
  }

  const_iterator cbegin() const {
    return impl_ptr_->cbegin();
  }

  const_iterator cend() const {
    return impl_ptr_->cend();
  }

  const_reverse_iterator rbegin() const {
    return impl_ptr_->rbegin();
  }

  const_reverse_iterator rend() const {
    return impl_ptr_->rend();
  }

  //! Returns reduced version of this word.
  Word freelyReduce() const {
    return freelyReduce(begin(), end());
  }

  //! Returns a copy of this word with reduced segment [begin, end).
  Word freelyReduce(const_iterator begin, const_iterator end) const;
  // freely reduces a segment of a word specified by [begin, end].
  // Most of the operations produce freely reduced words,
  // except operations insert, remove. The function takes linear time to performs
  // which is not efficient. So, either try to avoid uses of functions insert and remove
  // or try to bound the changes in variables beg and end

  //! Generates a pseudo randomly reduced word of the length wLen over the alphabet \f$X = \{x_1,\ldots,x_n\}\f$.
  static Word randomWord(int n, int wLen);

  //! Generates a pseudo randomly reduced word of a length in [wLenMin,wLenMax] and  over the alphabet \f$X = \{x_1,\ldots,x_n\}\f$.
  static Word randomWord(int n, int wLenMin, int wLenMax);

  std::list<int> toList() const {
    return impl_ptr_->toList();
  }

  const std::vector<int>& toVector() const {
    return impl_ptr_->toVector();
  }

    //! Get the length of the word.
  inline size_t length() const {
    return impl_ptr_->length();
  }

  //! Get the length of the word.
  inline size_t size() const {
    return impl_ptr_->size();
  }

  //! Checks if the word is empty.
  inline bool empty() const {
    return impl_ptr_->empty();
  }

  //! Multiply the word by a one-letter word defined by gen on the right. The result is being reduced.
  inline Word& push_back(int gen) {
    clone_();
    impl_ptr_->push_back(gen);
    return *this;
  }

  //! Multiply the word by a one-letter word defined by gen on the left. The result is being reduced.
  inline Word& push_front(int gen) {
    clone_();
    impl_ptr_->push_front(gen);
    return *this;
  }

  //! Multiply the word by a word on the right. The result is being reduced.
  Word& push_back(const Word& w);

  //! Multiply the word by a word on the left. The result is being reduced.
  Word& push_front(const Word& w);

  //! Remove the last symbol
  void pop_back() {
    clone_();
    impl_ptr_->pop_back();
  }

  //! Remove the first symbol
  void pop_front() {
    clone_();
    impl_ptr_->pop_front();
  }

  //! Determines the power of the word (as an element of a free monoid, not as an element of a free group).
  int getPower(Word& base) const {
    WordRep b;
    int exp;

    std::tie(b, exp) = impl_ptr_->root();

    base = Word(std::move(b));

    return exp;
  }


  //! Checks if the word contains a generator given by gen.
  bool doesContain(int gen) const {
    return impl_ptr_->contains(gen);
  }

  //! Shifts the word one position to the left, i.e., \f$cyclicLeftShift( x_1 x_2 \ldots x_{k-1} x_k ) = x_2 \ldots x_{k-1} x_k x_1 \f$.
  inline void cyclicLeftShift() {
    clone_();
    impl_ptr_->cyclicLeftShift();
  }

  //! Shifts the word one position to the right, i.e., \f$cyclicRightShift( x_1 x_2 \ldots x_{k-1} x_k ) = x_k x_1 x_2 \ldots x_{k-1} \f$.
  inline void cyclicRightShift() {
    clone_();
    impl_ptr_->cyclicRightShift();
  }

  //! Returns the cyclically reduced word.
  inline Word cyclicallyReduce() const {
    auto result = clone();

    result.impl_ptr_->cyclicallyReduce();

    return result;
  }

  //! Cyclically reduces the Word
  inline void cyclicallyReduceWord() {
    clone_();
    impl_ptr_->cyclicallyReduce();
  }

  //! Returns the cyclically reduced word and the corresponding conjugator.
  Word cyclicallyReduce(Word& conjugator) const {
    auto result = clone();

    auto c = result.impl_ptr_->cyclicallyReduce();
    conjugator = Word(std::move(c));

    return result;
  }

  //! Cyclically reduces the Word and returns the corresponding conjugator.
  void cyclicallyReduceWord(Word& conjugator) {
    clone_();

    auto c = impl_ptr_->cyclicallyReduce();
    conjugator = Word(std::move(c));
  }

  //! Inverts the word.
  inline Word inverse() const {
    return Word(impl_ptr_->inverse());
  }

  //! Cyclically permutes the word and return the result.
  /*! n>0 => left-shift, n<0 => rigth-shift permute. */
  inline Word cyclicallyPermute(int n) const {
    auto result = clone();
    result.impl_ptr_->cyclicallyPermute(n);
    return result;
  }

  //! Cyclically permutes the word.
  inline Word& cyclicallyPermuteWord(int n) {
    clone_();
    impl_ptr_->cyclicallyPermute(n);
    return *this;
  }

  //! Gets an initial segment [0, to) of the word.
  inline Word initialSegment(size_t to) const {
    auto result = clone();
    result.impl_ptr_->initialSegment(to);
    return result;
  }

  //! Gets a terminal segment [from, size) of the word.
  inline Word terminalSegment(size_t from) const {
    auto result = clone();
    result.impl_ptr_->terminalSegment(from);
    return result;
  }

  // Gets a segment [from, to) of the word.
  inline Word segment(size_t from, size_t to) const {
    auto result = clone();
    result.impl_ptr_->segment(from, to);
    return result;
  }

  Word subword(size_t from, size_t to) const {
    return Word(impl_ptr_->subword(from, to));
  }

  inline int exponentSum(int gen) const {
    return impl_ptr_->exponentSum(gen);
  }

  int isIn(int gen) const {
    return impl_ptr_->occurrences(gen);
  }

  Word power(int t) const;

  //! Insert a sequence of generators [b, e) into the word at a position pos
  template<class InputIterator>
  void insert(size_t pos, InputIterator b, InputIterator e) {
    clone_();
    impl_ptr_->insert(pos, b, e);
  }

  //! Insert a generator g into the word at a position pos
  void insert(size_t pos, int g);

  //! Replace a generator at a position pos by g
  void replace(size_t pos, int g);

  //! Replaces a subword of a word starting at a position it by a word [b, e).
  //! The length of the word does not increase if [b, e) is longer than the terminal segment of the word [it, end()).
  //! In that case terminal symbols of [b, e) are ignored.
  template<class InputIterator>
  void replace(size_t pos, InputIterator b, InputIterator e) {
    clone_();
    impl_ptr_->replace(pos, b, e);
  }

  //! Returns a word in which generators are replaced by words
  /*! Replace generators with the words (images) contained 
    in the vector \param images
    */
  Word replaceGenerators(const std::vector<Word>& images) const;

  //! Compute the minimal equivalent word.
  /*! permutableGenerators must be positive.
   */
  Word minimalEquivalentForm(const std::set<int>& permutableGenerators, bool inverses, bool cyclicPermutations) const;

  std::string toString() const;

private:
  friend ostream& operator<<(ostream& os, const Word& w) {
    InfiniteAlphabet::defaultAlphabet.printWord(os, w);
    return os;
    //return w.look( )->printOn( os );
  }

  friend istream& operator>>(istream& is, Word& w) {
    w = InfiniteAlphabet::defaultAlphabet.readWord(is);
    return is;
    //return w.look( )->printOn( os );
  }

  ostream& printOn(ostream& os) const {
    InfiniteAlphabet::defaultAlphabet.printWord(os, *this);
    return os;
    //return look( )->printOn( os );
  }

  //! Makes a deep copy of this word.
  Word clone() const {
    auto result = *this;
    result.clone_();
    return result;
  }

private:
  //! Clones the underlying word rep.
  void clone_() {
    WordRep copy = *impl_ptr_;
    impl_ptr_ = std::make_shared<WordRep>(std::move(copy));
  }

  std::shared_ptr<WordRep> impl_ptr_;
};

Word operator"" _w(const char* str, size_t);

Word abelianization(const Word& w);

//! For a word w returns a map m such that for i > 0
//! m[i] is the number of occurrences of x_{i}^{+/- 1} in w.
std::map<size_t, size_t> occurrences(const Word& w);

#endif // CRAG_WORD_H
