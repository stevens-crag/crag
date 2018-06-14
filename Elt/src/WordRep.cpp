// Copyright (C) 2005 Alexander Ushakov
// Contents: Implementation of class WordRep
//
// Principal Authors: Alexander Ushakov
//


#include "WordRep.h"

#include <cmath>
#include <sstream>

WordRep::WordRep(const std::list<int>& gens)
  : elements_(gens.begin(), gens.end()) {
  validate_();
  freelyReduce();
}

WordRep::WordRep(std::vector<int> gens)
  : elements_(std::move(gens)) {
  validate_();
  freelyReduce();
}


WordRep::WordRep(int g)
  : elements_({g}) {
  validate_();
}

WordRep::WordRep(std::initializer_list<int> gens)
  : elements_(gens) {
  validate_();
  freelyReduce();
}

std::string WordRep::toString() const {
  std::ostringstream os;

  os << *this;

  return os.str();
}

WordRep::iterator WordRep::begin() {
  return elements_.begin();
}

WordRep::const_iterator WordRep::begin() const {
  return elements_.begin();
}

WordRep::iterator WordRep::end() {
  return elements_.end();
}

WordRep::const_iterator WordRep::end() const {
  return elements_.end();
}

WordRep::const_iterator WordRep::cbegin() const {
  return elements_.cbegin();
}

WordRep::const_iterator WordRep::cend() const {
  return elements_.cend();
}

WordRep::const_reverse_iterator WordRep::rbegin() const {
  return elements_.rbegin();
}

WordRep::const_reverse_iterator WordRep::rend() const {
  return elements_.rend();
}

WordRep& WordRep::operator^=(int power) {
  const auto base = (power > 0) ? *this : this->inverse();

  this->clear();

  for (int i = 0; i < std::abs(power); ++i) {
    *this *= base;
  }

  return *this;
}

WordRep& WordRep::operator^=(const WordRep& conjugator) {
  auto result = conjugator.inverse();
  result *= *this;
  result *= conjugator;

  elements_.swap(result.elements_);

  return *this;
}

WordRep& WordRep::operator*=(const WordRep& other) {
  for (const auto g : other.elements_) {
    reduced_push_back_(g);
  }

  return *this;
}

bool WordRep::contains(int gen) const {
  validate_(gen);

  return std::any_of(elements_.begin(), elements_.end(), [gen](int el) {
    return std::abs(el) == std::abs(gen);
  });
}

void WordRep::cyclicLeftShift() {
  if (elements_.size() < 2) {
    return;
  }

  reduced_push_back_(elements_.front());
  elements_.erase(elements_.begin());
}

void WordRep::cyclicRightShift() {
  if (elements_.size() < 2) {
    return;
  }

  reduced_push_front_(elements_.back());
  elements_.pop_back();
}


//=========================================================


size_t WordRep::occurrences(int gen) const {
  if (gen == 0) {
    throw std::invalid_argument("Zero indices are not allowed.");
  }

  const int abs_gen = std::abs(gen);
  size_t result = 0;

  for (const auto el : elements_) {
    if (abs_gen == std::abs(el)) {
      ++result;
    }
  }

  return result;
}

int WordRep::exponentSum(int gen) const {
  if (gen == 0) {
    throw std::invalid_argument("Zero indices are not allowed.");
  }

  const int abs_gen = std::abs(gen);
  int result = 0;

  for (const auto el : elements_) {
    const int degree = (el > 0) ? 1 : -1;

    if (abs_gen == std::abs(el)) {
      result += degree;
    }
  }

  return result;
}

WordRep WordRep::inverse() const {
  storage_t elements;
  elements.reserve(size());

  for (auto it = elements_.rbegin(); it != elements_.rend(); ++it) {
    elements.push_back(-*it);
  }

  return WordRep(std::move(elements));
}

bool WordRep::operator==(const WordRep& other) const {
  return elements_ == other.elements_;
}

bool WordRep::operator!=(const WordRep& other) const {
  return !(*this == other);
}

bool WordRep::operator<(const WordRep& other) const {
  const auto this_size = size();
  const auto other_size = other.size();

  if (this_size != other_size) {
    return this_size < other_size;
  }

  auto this_it = elements_.begin();
  auto other_it = other.elements_.begin();

  for (size_t t = 0; t < this_size; ++t) {
    if (*this_it < *other_it) {
      return true;
    }

    if (*(this_it++) > *(other_it++)) {
      return false;
    }
  }

  // words are equal
  return false;
}

bool WordRep::operator>(const WordRep& other) const {
  return other < *this;
}

void WordRep::insert(iterator position, int g) {
  validate_({g});

  elements_.insert(position, g);
}

void WordRep::insert(size_t position, int g) {
  position = std::min(position, size());

  auto it = elements_.begin();
  std::advance(it, position);

  insert(it, g);
}

void WordRep::replace(iterator position, int g) {
  validate_({g});

  *position = g;
}

void WordRep::replace(size_t position, int g) {
  if (position >= size()) {
    throw std::invalid_argument("Bad position.");
  }

  auto it = elements_.begin();
  std::advance(it, position);

  replace(it, g);
}

void WordRep::cyclicallyPermute(int n) {
  if (size() < 2) {
    return;
  }

  int sz = static_cast<int>(size());

  // make n positive
  n = ((n % sz) + sz) % sz;

  if (n == 0) {
    return;
  }

  auto middle = elements_.begin();
  std::advance(middle, n);

  std::rotate(elements_.begin(), middle, elements_.end());

  freelyReduce();
}

void WordRep::segment(size_t from, size_t to) {
  if (from > to) {
    throw std::invalid_argument("Bad segment.");
  }

  initialSegment(to);
  terminalSegment(from);
}

WordRep WordRep::subword(size_t from, size_t to) const {
  if (from > to) {
    throw std::invalid_argument("Bad subword.");
  }

  from = std::min(from, size());
  to = std::min(to, size());

  auto it_begin = elements_.begin();
  std::advance(it_begin, from);

  auto it_end = elements_.begin();
  std::advance(it_end, to);

  return WordRep(it_begin, it_end);
}

void WordRep::initialSegment(size_t to) {
  to = std::min(to, size());

  elements_.resize(to);
}

void WordRep::terminalSegment(size_t from) {
  if (from > size()) {
    from = size();
  }

  auto it = elements_.begin();
  std::advance(it, from);

  elements_.erase(elements_.begin(), it);
}

WordRep WordRep::cyclicallyReduce() {
  WordRep conjugator;

  // TODO: rewrite using iterators to avoid erase
  while (elements_.size() > 1) {
    const auto b = elements_.front();
    const auto e = elements_.back();

    if (e + b != 0) {
      break;
    }

    elements_.pop_back();
    elements_.erase(elements_.begin());

    conjugator.elements_.insert(conjugator.elements_.begin(), e);
  }

  return conjugator;
}

std::pair<WordRep, int> WordRep::root() const {
  auto len = size();

  if (len <= 1) {
    return std::make_pair(*this, 1);
  }

  int up = std::sqrt(static_cast<double>(len)) + 1;

  // Decompose len into a product of primes
  std::list< std::pair< int , int > > primes;
  std::list< std::pair< int , int > > primes_to_use;

  for( int i=2 ; i<up ; ++i ) {
    bool prime = true;
    auto p_it = primes.begin( );
    for( ; p_it!=primes.end( ) ; ++p_it ) {
      int a = (*p_it).first;
      int b = (*p_it).second;
      (*p_it).second = b+1==a ? (prime=false) : b+1;
    }
    if( prime ) {
      primes.push_back( std::pair< int , int >(i,0) );
      int count = 0;
      while( len%i==0 ) {
        ++count;
        len /= i;
      }
      if( count ) {
        primes_to_use.push_back( std::pair<int,int>(i,count) );
        up = std::sqrt(static_cast<double>(len)) + 1;
      }
    }
  }
  if( len>1 )
    primes_to_use.push_back( std::pair<int,int>(len,1) );

  /*{
	list< pair< int , int > >::iterator p_it = primes_to_use.begin( );
	for( ; p_it!=primes_to_use.end( ) ; ++p_it ) {
		int a = (*p_it).first;
		int b = (*p_it).second;
		cout << " " << a << "," << b << endl;
	}
	}*/

  len = size();

  // Construct a vector of generators (to have indexed access)
  std::vector< int > wrd_vct(len);
  auto g_it = cbegin( );
  for( int j=0 ; g_it!=cend() ; ++j, ++g_it )
    wrd_vct[j] = *g_it;

  // Find the power
  auto p_it = primes_to_use.begin( );
  for( ; p_it!=primes_to_use.end( ) ; ++p_it ) {
    int parts = (*p_it).first;
    int count = (*p_it).second;
    int shrt_len = len/parts;

    bool progress = true;
    for( int t=0 ; t<count && progress ; ++t ) {
      for( int offset=0 ; offset<shrt_len  && progress ; ++offset ) {
        int e = wrd_vct[offset];
        for( int j=offset ; j<len && progress ; j+=shrt_len )
          if( wrd_vct[j]!=e )
            progress = false;
      }
      if( progress )
        shrt_len = (len = shrt_len) / parts;
    }
  }

  auto base = *this;
  base.initialSegment( len );
  return std::make_pair(base, size() / len);
};

void WordRep::freelyReduce() {
  for (auto it = elements_.begin(); it != elements_.end();) {
    auto next = it;
    ++next;

    if (next == elements_.end()) {
      return;
    }

    if (*it + *next != 0) {
      it = next;
      continue;
    }

    // reduction
    const auto is_beginning = (it == elements_.begin());

    it = elements_.erase(elements_.erase(it));

    if (!is_beginning) {
      --it;
    }
  }
}

void WordRep::freelyReduce(iterator begin, iterator end) {
  // cut suffix
  storage_t suffix(end, elements_.end());
  elements_.erase(end, elements_.end());

  // cut prefix
  storage_t prefix(elements_.begin(), begin);
  elements_.erase(elements_.begin(), begin);


  // reduce the remaining part
  freelyReduce();

  // put prefix and suffix back
  elements_.insert(elements_.begin(), prefix.begin(), prefix.end());
  elements_.insert(elements_.end(), suffix.begin(), suffix.end());
}

void WordRep::reduced_push_back_(int g) {
  validate_({g});

  if (elements_.empty()) {
    elements_.push_back(g);
    return;
  }

  const int last_el = elements_.back();

  if (last_el + g == 0) {
    elements_.pop_back();
  } else {
    elements_.push_back(g);
  }
}

void WordRep::reduced_push_front_(int g) {
  validate_({g});

  if (elements_.empty()) {
    elements_.push_back(g);
    return;
  }

  const int first_el = elements_.front();

  if (first_el + g == 0) {
    elements_.erase(elements_.begin());
  } else {
    elements_.insert(elements_.begin(), g);
  }
}

std::ostream& operator<<(std::ostream& out, const WordRep& w) {
  if (w.empty()) {
    return out << "1";
  }

  const std::string letter("x");

  int last_index = w.front();
  int last_degree = 1;

  const auto printLastEl = [&] () {
    out << letter << std::abs(last_index);

    if (last_index < 0 || last_degree != 1) {
      out << "^" << ((last_index > 0) ? last_degree : -last_degree);
    }
  };

  for (auto it = ++w.begin(); it != w.end(); ++it) {
    const int index = *it;

    if (last_index != index) {
      printLastEl();

      if (last_degree != 0) {
        out << " ";
      }

      last_index = index;
      last_degree = 1;
    } else {
      ++last_degree;
    }
  }

  printLastEl();

  return out;
}
