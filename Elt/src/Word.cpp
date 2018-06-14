// Copyright (C) 2005 Alexander Ushakov
// Contents: Implementation of class Word
//
// Principal Authors: Alexander Ushakov
//

#include "Word.h"

#include <map>

#include "RanlibCPP.h"

Word& Word::push_back(const Word& w) {
  clone_();

  for (const auto g : w) {
    impl_ptr_->push_back(g);
  }

  return *this;
}

Word& Word::push_front(const Word& w) {
  clone_();

  auto it = w.end();

  while (it != w.begin()) {
    impl_ptr_->push_front(*--it);
  }

  return *this;
}

Word Word::power(int t) const {
  if (t == 0) {
    return Word();
  }

  auto result = clone();
  *result.impl_ptr_ ^= t;
  return result;
}

Word Word::randomWord(int gens, int wLen) {
  if (wLen == 0) return Word();

  int old = 0;
  std::list<int> result;
  for (int i = 0; i < wLen; ++i) {
    int div = i == 0 ? 2 * gens : 2 * gens - 1;
    // int g = ::rand()%div-gens;
    int g = RandLib::ur.irand(0, div - 1) - gens;
    g = g >= 0 ? g + 1 : g;
    if (g + old == 0)
      g = gens;
    result.push_back(old = g);
  }

  return Word(std::move(result));
}

Word Word::randomWord(int gens, int wLenMin, int wLenMax) {

  int wLen = RandLib::ur.irand(wLenMin, wLenMax);
  return randomWord(gens, wLen);
}

Word Word::freelyReduce(const_iterator begin, const_iterator end) const {
  const auto from = std::distance(this->begin(), begin);
  const auto to = std::distance(this->begin(), end);

  auto result = clone();

  auto new_begin = result.impl_ptr_->begin();
  std::advance(new_begin, from);

  auto new_end = result.impl_ptr_->begin();
  std::advance(new_end, to);

  result.impl_ptr_->freelyReduce(new_begin, new_end);

  return result;
}

void Word::insert(size_t pos, int g) {
  clone_();
  impl_ptr_->insert(pos, g);
}

void Word::replace(size_t pos, int g) {
  clone_();
  impl_ptr_->replace(pos, g);
}

Word Word::replaceGenerators(const std::vector<Word>& images) const {
  Word result;

  for (const auto g : *this) {
    const auto g_abs = std::abs(g);

    if (g_abs > images.size()) {
      throw std::invalid_argument("Word::replaceGenerators() : image vector index overflow.");
    }

    result *= (g > 0 ? images[abs(g) - 1] : images[abs(g) - 1].inverse());
  }

  return result;
}

Word Word::minimalEquivalentForm(const std::set<int>& permutableGenerators, bool inverses, bool cyclicPermutations) const {
  Word w;
  int P = 1;
  if (cyclicPermutations) {
    w = this->cyclicallyReduce();
    P = w.length();
  } else {
    w = *this;
  }
  Word result = w;
  int I = inverses ? 2 : 1;

  // cout << "   >>>   " << cur << endl;

  for (int i = 0; i < I; ++i) {
    Word cur = i == 0 ? w : -w;
    for (int p = 0; p < P; ++p, cur.cyclicLeftShift()) {

      // prepare a permutation
      std::map<int, int> permutation;
      for (std::set<int>::const_iterator g_it = permutableGenerators.begin(); g_it != permutableGenerators.end(); ++g_it)
        permutation[*g_it] = permutation[-*g_it] = 0;

      // slide along the word and replace permutable generators to guarantee minimality
      Word rep;
      std::set<int>::const_iterator counter_it = permutableGenerators.end();
      for (auto w_it = cur.cbegin(); w_it != cur.cend(); ++w_it) {
        int g = *w_it;
        // check if g is permutable
        std::map<int, int>::iterator p_it = permutation.find(g);
        if (p_it == permutation.end()) { // non-permutable
          rep.push_back(g);
        } else { // permutable
          if ((*p_it).second == 0) {
            (*p_it).second = -*(--counter_it);
            permutation[-g] = *counter_it;
          }
          rep.push_back((*p_it).second);
        }
      }
      // cout << rep << endl;
      if (rep < result)
        result = rep;

      //for( set< int >::const_iterator g_it = permutableGenerators.begin( ) ; g_it!=permutableGenerators.end( ) ; ++g_it )
      // cout << "  " << permutation[*g_it] << " , " << permutation[-*g_it] << endl;

    }
  }

  // cout << "   <<<   " << result << endl;

  return result;
}

std::string Word::toString() const {
  std::ostringstream os;

  os << *this;

  return os.str();
}

Word operator"" _w(const char* str, size_t) {
  Word w;
  std::istringstream s(str);

  s >> w;

  return w;
}

Word abelianization(const Word& w) {
  std::map<int, int> m;

  for (const auto l : w) {
    if (l > 0) {
      m[l] += 1;
    } else {
      m[-l] -= 1;
    }
  }

  std::list<int> result;

  for (const auto& p : m) {
    const auto gen = p.first;
    const auto exp = p.second;

    if (exp == 0) {
      continue;
    }

    result.insert(result.end(), std::abs(exp), (exp > 0) ? gen : -gen);
  }

  return Word(std::move(result));
}

std::map<size_t, size_t> occurrences(const Word& w) {
  std::map<size_t, size_t> result;

  for (const auto gen : w) {
    const auto abs_gen = std::abs(gen);

    if (result.count(abs_gen) == 0) {
      result[abs_gen] = 0;
    }

    result[abs_gen] += 1;
  }

  return result;
}
