// Copyright (C) 2018 Alexander Ushakov
// Contents: Implementation of class Permutation
//
// Principal Authors: Alexander Ushakov
//

#include "Permutation.h"

#include <cassert>
#include <limits>
#include <list>
#include <set>
#include <sstream>

#include "RanlibCPP.h"

Permutation::Permutation(size_t l)
    : values_(l) {
  for (size_t i = 0; i < l; ++i) {
    values_[i] = i;
  }
}

Permutation::Permutation(std::vector<int> values)
    : values_(std::move(values)) {
  validate_();
}

std::string Permutation::toString() const {
  std::ostringstream os;

  os << *this;

  return os.str();
}

Permutation Permutation::power(int p) const {
  Permutation base = (p < 0) ? -*this : *this;

  Permutation result(size());

  auto ap = std::abs(p);

  while (ap) {
    if (ap % 2 == 1) {
      result *= base;
    }

    ap = ap >> 1;

    base = base * base;
  }

  return result;
}

Permutation& Permutation::left_mult_by_cycle(const std::vector<int>& cycle) {
  for (size_t i = 1; i < cycle.size(); ++i)
    std::swap(values_[cycle[i - 1]], values_[cycle[i]]);

  return *this;
}

void Permutation::lr_multiply_by_cycles(
    Permutation& P, Permutation& I, const std::vector<int>& M1, const std::vector<int>& M2) {
  for (size_t i = 1; i < M1.size(); ++i) {
    int a = M1[i - 1];
    int b = M1[i];
    int c = P.values_[a];
    int d = P.values_[b];
    std::swap(P.values_[a], P.values_[b]);
    std::swap(I.values_[c], I.values_[d]);
  }

  for (int i = M2.size() - 1; i > 0; --i) {
    int a = M2[i - 1];
    int b = M2[i];
    int c = I.values_[a];
    int d = I.values_[b];
    std::swap(I.values_[a], I.values_[b]);
    std::swap(P.values_[c], P.values_[d]);
  }
}

Permutation Permutation::computeConjugacyClassRepresentative(Permutation& conj) const {
  int i;
  int size = values_.size();
  Permutation p1(*this);
  conj = Permutation(size);

  // 1. arrange cycles
  std::set<std::pair<int, int>> cycles;
  for (i = 0; i < size; ++i) {
    if (p1.values_[i] != i) {
      int r = i;
      do {
        if (p1.values_[r] != r + 1) {
          Permutation c(size);
          c.change(r + 1, p1.values_[r]);

          p1 = c * p1 * c;
          conj *= c;
        }
        r = p1.values_[r];
      } while (p1.values_[r] != i);
      cycles.insert(std::pair<int, int>(r + 1 - i, i));
      i = r;
    } else
      cycles.insert(std::pair<int, int>(1, i));
  }

  auto c_it = cycles.begin();

  int pos = 0;
  while (cycles.size()) {
    c_it = cycles.begin();
    int len = (*c_it).first;
    int pt = (*c_it).second;
    cycles.erase(c_it);

    if (pt != pos) {
      // 2.1. prepare shift conjugators
      Permutation c1(size);
      Permutation c2(size);
      Permutation c3(size);

      for (i = pos; i < pt + len; ++i)
        c1.values_[i] = pt + len - 1 - (i - pos);

      for (i = pos; i < pos + len; ++i)
        c2.values_[i] = pos + len - 1 - (i - pos);

      for (i = pos + len; i < pt + len; ++i)
        c3.values_[i] = pt + len - 1 - (i - pos - len);

      // 2.2. shift cycle to a position pos
      p1 = c1 * p1 * c1.inverse();
      p1 = c2 * p1 * c2.inverse();
      p1 = c3 * p1 * c3.inverse();
      conj *= c1.inverse();
      conj *= c2.inverse();
      conj *= c3.inverse();

      std::set<std::pair<int, int>> _cycles;
      c_it = cycles.begin();
      for (; c_it != cycles.end(); ++c_it) {
        if ((*c_it).second < pt)
          _cycles.insert(std::pair<int, int>((*c_it).first, (*c_it).second + len));
        else
          _cycles.insert(std::pair<int, int>((*c_it).first, (*c_it).second));
      }
      cycles = _cycles;
    }
    pos += len;
  }

  return p1;
}

Permutation Permutation::computeConjugator(const Permutation& p) const {
  int size = values_.size();
  assert(size == p.values_.size());
  //  if( size!=p.theValue.size( ) ) Do something with this!!!
  //    return false;

  Permutation result(size);

  Permutation conj1(size);
  Permutation p1 = computeConjugacyClassRepresentative(conj1);

  Permutation conj2(size);
  Permutation p2 = p.computeConjugacyClassRepresentative(conj2);

  // cout << "_p1 = " << p1 << endl;
  // cout << "_p2 = " << p2 << endl;

  result = conj1 * conj2.inverse();
  return result;
}

/*
Permutation
Permutation::computeConjugator( const Permutation& p ) const
{
  int size = theValue.size( );
  if( size!=p.theValue.size( ) )
    return false;

  Permutation result( size );

  vector< bool > h( size , false );
  for( int i=size-1 ; i>=0 ; --i ) {
    int r = i;
    int s = i;
    while( !h[r] ) {
      h[r] = true;
      r = theValue[r];
      s = p.theValue[s];
      if( r!=s )
        result.theValue[r] = s;
    }
  }

  // cout << " ---- " << result << endl;

  vector< bool > k( size , false );
  for( int i=size-1 ; i>=0 ; --i ) {
    if( !k[i] ) {
      k[i] = true;
      int r = i;
      while( !k[result.theValue[r]] && result.theValue[r]!=r ) {
        r = result.theValue[r];
        k[r] = true;
      }

      result.theValue[r] = i;
    }
  }

  // cout << " ---- " << result << endl;

  return result;
}
*/

/*
bool
Permutation::computeConjugator( const Permutation& p , Permutation& res ) const
{
  Permutation _p( *this );

  int size = _p.theValue.size( );
  if( size!=p.theValue.size( ) )
    return false;

  res = Permutation( size );

  cout << "------------------" << endl;
  cout << _p << endl;
  cout << p << endl;

  for( int i=0 ; i<size ; ++i ) {

    int y = p.theValue[i];
    Permutation inv = _p.inverse( );
    int k = inv[y];
    if( k==i )
      continue;

    cout << "i = " << i << endl;
    cout << "k = " << k << endl;
    Permutation c( size );
    c.change( i , k );

    _p = c * _p * c;

    if( _p.theValue[i]!=y )
      return false;
    cout << "------------------" << endl;
    cout << _p << endl;
    cout << p << endl;

    res *= c;
  }

  return true;
}
*/

std::vector<int> Permutation::getWordPresentation() const {
  int j;
  std::vector<int> result;

  /*
  // first implementation of this function
  // much worse than the second one
  vector<int> val( theValue );
  for( int i=0 ; i<val.size() ; ++i ) {
    if( val[i]!=i ) {
      int t = val[i];
      int s = i;

      for( int j=t-1 ; j>s ; --j )
        result.push_back( -j-1 );
      for( int j=s ; j<t ; ++j )
        result.push_back( j+1 );

      swap( val[i] , val[val[i]] );
      i--;
    }
  }
  */


  // second implementation
  std::vector<int> val(values_);
  Permutation inv = inverse();
  std::vector<int> inv_val(inv.values_);
  for (int i = 0; i < values_.size(); ++i) {
    if (val[i] != i) {
      int r = val[i];
      int t = inv_val[i];
      int s = i;

      for (j = t - 1; j > s; --j)
        result.push_back(-j - 1);
      for (j = s; j < t; ++j)
        result.push_back(j + 1);

      std::swap(val[s], val[t]);
      std::swap(inv_val[s], inv_val[r]);
    }
  }

  return result;
}


std::vector<int> Permutation::geodesic() const {
  std::vector<int> result;

  Permutation cur(size());
  Permutation inv(size());

  for (int i = 0; i < size(); ++i) {
    int pos = inv.values_[values_[i]];

    for (int j = pos - 1; j >= i; --j) {
      result.push_back(j);
      inv.change(cur.values_[j], cur.values_[j + 1]);
      cur.change(j, j + 1);
    }
  }

  for (size_t i = 0; i < result.size() / 2; ++i) {
    std::swap(result[i], result[result.size() - i - 1]);
  }

  return result;
}

std::vector<int> Permutation::geodesicWord() const {
  auto vec = geodesic();

  for (size_t i = 0; i < vec.size(); ++i) {
    ++vec[i];
  }

  return vec;
}

bool Permutation::isTrivial() const {
  for (size_t i = 0; i < values_.size(); ++i) {
    if (values_[i] != i) {
      return false;
    }
  }

  return true;
}

bool Permutation::operator==(const Permutation& p) const {
  return values_ == p.values_;
}

bool Permutation::operator!=(const Permutation& p) const {
  return values_ != p.values_;
}

bool Permutation::operator<(const Permutation& p) const {
  if (size() != p.size()) {
    return size() < p.size();
  }

  for (size_t t = 0; t < values_.size(); ++t) {
    if (values_[t] < p.values_[t]) {
      return true;
    }

    if (values_[t] > p.values_[t]) {
      return false;
    }
  }

  return false;
}

Permutation Permutation::operator*(const Permutation& p) const {
  return (Permutation(*this) *= p);
}

Permutation& Permutation::operator*=(const Permutation& other) {
  if (size() < other.size()) {
    for (size_t i = size(); i < other.size(); ++i) {
      values_.push_back(i);
    }
  }

  for (size_t t = 0; t < size(); ++t) {
    if (values_[t] < other.size()) {
      values_[t] = other[values_[t]];
    }
  }

  return *this;
}

Permutation Permutation::inverse() const {
  std::vector<int> result(size());

  for (size_t t = 0; t < size(); ++t) {
    result[values_[t]] = t;
  }

  return Permutation(std::move(result));
}

size_t Permutation::difference(const Permutation& p) const {
  if (size() != p.size()) {
    return std::numeric_limits<size_t>::max();
  }

  size_t result = 0;

  for (size_t i = 0; i < size(); ++i) {
    if (values_[i] != p[i]) {
      ++result;
    }
  }

  return result;
}

Permutation Permutation::random(size_t n) {
  Permutation res(n);

  for (size_t t = 0; t + 1 < n; ++t) {
    auto pos = RandLib::ur.irand(0, n - t - 1);
    std::swap(res.values_[t], res.values_[t + pos]);
  }

  return res;
}

std::ostream& operator<<(std::ostream& os, const Permutation& p) {
  const auto len = p.size();

  os << "{";

  for (size_t t = 0; t < len; ++t) {
    if (t > 0) {
      os << ", ";
    }

    os << p[t];
  }

  os << "}";

  return os;
}

Permutation Permutation::join2(const Permutation& p) const {
  int i;
  int l1 = values_.size();
  int l2 = p.values_.size();
  if (l1 != l2) {
    std::cerr << "Check dimensions in join2 operation" << std::endl;
    exit(1);
  }


  /*
  Permutation omega( l1 );
  for( int i=0 ; i<l1 ; ++i )
    omega.theValue[i] = (i+l1-1)%l1;

  Permutation p1 = inverse( )*omega;
  Permutation p2 = p.inverse( )*omega;
  Permutation p3 = p1.meet( p2 );

  return p3.inverse( )*omega;
  */

  std::vector<int> cycleN(l1, -1);
  std::vector<int> foundPts(l1, -1);
  for (i = l1 - 1; i >= 0; --i) {
    if (cycleN[i] != -1)
      continue;

    std::list<int> toCheck;
    toCheck.push_back(i);
    foundPts[i] = i;
    while (toCheck.begin() != toCheck.end()) {
      int cur = *toCheck.begin();
      toCheck.pop_front();
      cycleN[cur] = i;

      int next1 = values_[cur];
      int next2 = p.values_[cur];

      if (foundPts[next1] != i) {
        foundPts[next1] = i;
        toCheck.push_back(next1);
      }
      if (foundPts[next2] != i) {
        foundPts[next2] = i;
        toCheck.push_back(next2);
      }
    }
  }

  std::vector<int> perm(l1, -1);

  std::vector<int> firstPt(l1, -1);
  std::vector<int> prevPt(l1, -1);
  for (i = l1 - 1; i >= 0; --i) {
    int c = cycleN[i];
    perm[i] = i;
    if (firstPt[c] == -1)
      firstPt[c] = i;
    else
      std::swap(perm[i], perm[prevPt[c]]);
    prevPt[c] = i;
  }

  return Permutation(perm);
}


Permutation Permutation::meet2(const Permutation& p) const {
  int i;
  int l1 = values_.size();
  int l2 = p.values_.size();
  if (l1 != l2) {
    std::cerr << "Check dimensions in meet2 operation" << std::endl;
    exit(1);
  }

  // 1. extract cycles from the second permutation
  std::vector<int> cycle2N(l1, 0);
  for (i = 0; i < l1; ++i) {
    if (p.values_[i] != i && cycle2N[i] == 0) {
      cycle2N[i] = i + 1;
      for (int t = p.values_[i]; t != i; t = p.values_[t])
        cycle2N[t] = i + 1;
    }
    // cout << cycle2N[i] << ",";
  }
  // cout << endl;

  // 2. extract cycles from the first permutation
  //    and compute the result
  Permutation result(l1);
  std::vector<int> cycle1N(l1, 0);

  for (i = 0; i < l1; ++i) {
    if (values_[i] != i && cycle1N[i] == 0) {
      cycle1N[i] = i + 1;
      std::vector<std::pair<int, int>> pairs;
      if (cycle2N[i])
        pairs.push_back(std::pair<int, int>(cycle2N[i], i));
      for (int t = values_[i]; t != i; t = values_[t]) {
        cycle1N[t] = i + 1;
        if (cycle2N[t])
          pairs.push_back(std::pair<int, int>(cycle2N[t], t));
      }

      if (pairs.size() < 1)
        continue;

      // cout << "> " << result << endl;
      preparePairs_(l1 + 1, result, pairs);
      // cout << "< " << result << endl;
    }
    // cout << cycle1N[i] << ",";
  }
  // cout << endl;

  return result;
}


void Permutation::preparePairs_(int N, Permutation& P, std::vector<std::pair<int, int>>& pairs) const {
  std::set<std::pair<int, int>> pairs1;
  for (size_t i = 0; i < pairs.size(); ++i)
    pairs1.insert(pairs[i]);

  std::set<std::pair<int, int>>::reverse_iterator it = pairs1.rbegin();
  int len = 1;
  // int beg = (*it).second;
  std::pair<int, int> prev_pair = (*it);
  it++;
  for (; it != pairs1.rend(); ++it) {
    // cout << "(" << (*it).first << "," << (*it).second << ")" << endl;

    if (prev_pair.first == (*it).first) {
      P.change(prev_pair.second, (*it).second);
      // cout << prev_pair.second << "," << (*it).second << endl;
    } else {
      len = 0;
      // beg = (*it).second;
    }
    len++;
    prev_pair = (*it);
  }

  /*
  vector< int > counts( N , 0 );
  for( int i=0 ; i<pairs.size( ) ; ++i )
    counts[pairs[i].first]++;

  for( int i=1 ; i<N ; ++i )
    counts[i] += counts[i-1];

  vector< pair<int,int> > pairs1( pairs.size( ) );
  vector< int > cur_num( N , 0 );
  for( int i=0 ; i<pairs.size( ) ; ++i ) {
    pairs1[--counts[pairs[i].first]] = pairs[i];
    // counts[triples[i].c1]--;
  }

  for( int i=1 ; i<N ; ++i ) {
    if( count[i]-count[i-1]>1 ) {
      vector< int > counts( N , 0 );

    }
  }
  pairs = pairs1;
  */
}

Permutation Permutation::getHalfTwistPermutation(size_t n) {
  std::vector<int> result(n);

  for (size_t i = 0; i < n; ++i) {
    result[i] = n - i - 1;
  }

  return Permutation(result);
}

Permutation Permutation::RightGCD(const Permutation& p) const {
  const auto this_size = size();

  if (this_size != p.size()) {
    throw std::invalid_argument("Cannot compute RightGCD of permutation of different sizes.");
  }

  int* l_ind_a = new int[this_size];
  int* l_ind_b = new int[this_size];
  int* r_ind_a = new int[this_size];
  int* r_ind_b = new int[this_size];

  Permutation result(this_size);
  _sub_meet(p, inverse(), p.inverse(), result, l_ind_a, l_ind_b, r_ind_a, r_ind_b, 0, this_size);

  delete[] l_ind_a;
  delete[] l_ind_b;
  delete[] r_ind_a;
  delete[] r_ind_b;

  return result;
}

Permutation Permutation::RightLCM(const Permutation& p) const {
  Permutation Delta = getHalfTwistPermutation(size());
  return (*this * Delta).RightGCD(p * Delta) * Delta;
}

Permutation Permutation::LeftGCD(const Permutation& p) const {
  return -((-*this).RightGCD(-p));
}

Permutation Permutation::LeftLCM(const Permutation& p) const {
  Permutation Delta = getHalfTwistPermutation(size());
  Permutation P3 = Delta * inverse();
  Permutation P4 = Delta * p.inverse();
  return (Delta * P3.RightGCD(P4)).inverse();
}

void Permutation::validate_() const {
  std::set<int> values(values_.begin(), values_.end());

  if (values.size() != size()) {
    throw std::invalid_argument("Permutation contains duplicates.");
  }

  if (values.empty()) {
    return;
  }

  if ((*values.begin() != 0) || *values.rbegin() + 1 != size()) {
    throw std::invalid_argument("Permutation must be indexed from 0 to n-1.");
  }
}

void Permutation::_sub_meet(
    const Permutation& p,
    const Permutation& ip1,
    const Permutation& ip2,
    Permutation& cur,
    int* l_ind_a,
    int* l_ind_b,
    int* r_ind_a,
    int* r_ind_b,
    int beg,
    int end) const {
  if (end - beg == 1) {
    return;
  }

  const int middle = beg + (end - beg) / 2;

  // I. reorder left and right parts of permutation according to meet operation
  _sub_meet(p, ip1, ip2, cur, l_ind_a, l_ind_b, r_ind_a, r_ind_b, beg, middle);
  _sub_meet(p, ip1, ip2, cur, l_ind_a, l_ind_b, r_ind_a, r_ind_b, middle, end);


  // II. merge left and right parts

  // 1. find left indices
  for (int i = middle - 1; i >= beg; --i) {
    int a = ip1.values_[cur.values_[i]];
    int b = ip2.values_[cur.values_[i]];

    if (i == middle - 1) {
      l_ind_a[i] = a;
      l_ind_b[i] = b;
    } else {
      l_ind_a[i] = a < l_ind_a[i + 1] ? a : l_ind_a[i + 1];
      l_ind_b[i] = b < l_ind_b[i + 1] ? b : l_ind_b[i + 1];
    }
  }

  // 2. find right indices
  for (int i = middle; i < end; ++i) {
    int a = ip1.values_[cur.values_[i]];
    int b = ip2.values_[cur.values_[i]];

    if (i == middle) {
      r_ind_a[i] = a;
      r_ind_b[i] = b;
    } else {
      r_ind_a[i] = a < r_ind_a[i - 1] ? r_ind_a[i - 1] : a;
      r_ind_b[i] = b < r_ind_b[i - 1] ? r_ind_b[i - 1] : b;
    }
  }

  // 3. merge lists
  int i1 = beg;
  int i2 = middle;
  int* new_sublist = new int[end - beg];

  for (int i = 0; i < end - beg; ++i) {
    if (i1 == middle) {
      new_sublist[i] = cur.values_[i2++];
      continue;
    }

    if (i2 == end) {
      new_sublist[i] = cur.values_[i1++];
      continue;
    }

    if (l_ind_a[i1] > r_ind_a[i2] && l_ind_b[i1] > r_ind_b[i2]) {
      new_sublist[i] = cur.values_[i2++];
    } else {
      new_sublist[i] = cur.values_[i1++];
    }
  }

  for (int i = 0; i < end - beg; ++i) {
    cur.values_[beg + i] = new_sublist[i];
  }

  delete[] new_sublist;
}

Permutation Permutation::tinyFlip(int sh) const {
  int len = values_.size();
  Permutation result(len);

  if (sh < 0)
    sh = sh - (sh / len - 1) * len;

  for (int t = 0; t < len; ++t)
    result.values_[(t + sh) % len] = (values_[t] + sh) % len;

  return result;
}

bool Permutation::mixable(const Permutation& p1, const Permutation& p2) {
  int l1 = p1.values_.size();
  int l2 = p2.values_.size();
  if (l1 != l2) {
    std::cerr << "Check dimensions in mixable operation" << std::endl;
    exit(1);
  }

  Permutation ip1 = p1.inverse();
  for (int i = 1; i < l1; ++i)
    if (ip1.values_[i - 1] > ip1.values_[i] && p2.values_[i - 1] < p2.values_[i])
      return false;

  return true;
}

Permutation Permutation::flip() const {
  int sz = size();
  Permutation result(sz);

  for (int i = 0; i < sz; ++i)
    result.values_[sz - i - 1] = sz - values_[i] - 1;

  return result;
}

size_t Permutation::length() const {
  return geodesic().size();
}

Permutation Permutation::increaseSize(size_t n) const {
  if (n <= size()) {
    return *this;
  }

  auto result = values_;
  result.resize(n);

  for (size_t i = size(); i < n; ++i) {
    result[i] = i;
  }

  return Permutation(result);
}

Permutation Permutation::getCyclePermutation(size_t n) {
  std::vector<int> result(n);

  for (size_t i = 0; i + 1 < n; ++i) {
    result[i] = i + 1;
  }

  if (n > 0) {
    result[n - 1] = 0;
  }

  return Permutation(result);
}

std::vector<std::vector<int>> toCycles(const Permutation& p) {
  std::vector<std::vector<int>> cycles;

  std::vector<int> items = p.getVector();

  const int visited = -1;

  for (int cycle_start = 0; cycle_start < items.size(); ++cycle_start) {
    if (items[cycle_start] == visited) {
      continue;
    }

    // check cycle
    auto current = items[cycle_start];
    std::vector<int> cycle = {cycle_start};

    while (true) {
      // cycle finished
      if (current == cycle_start) {
        items[current] = visited;
        break;
      }

      cycle.push_back(current);

      int next = items[current];
      items[current] = visited;
      current = next;
    }

    if (cycle.size() > 1) {
      cycles.push_back(std::move(cycle));
    }
  }

  return cycles;
}
