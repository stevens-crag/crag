#pragma once

#ifndef CRAG_WALNUT_ATTACK_H
#define CRAG_WALNUT_ATTACK_H

#include <boost/optional.hpp>
#include <fstream>
#include <future>

#include "LinkedBraidStructure.h"
#include "fast_conjugacy_check.h"
#include "fast_identity_check.h"
#include "parallel.h"
#include "walnut.h"

namespace crag {
namespace walnut {

using FF = finitefield::ZZ<199>;
typedef crag::braidgroup::BraidHasher<FF>::braid_hash_t braid_hash_t;


// TODO: move to braid group file
int braidAbelianization(const Word& w) {
  int result = 0;

  for (const auto g : w) {
    result += (g > 0) ? 1 : -1;
  }

  return result;
}

// Big positive numbers must be interpreted as long cancellation, w1 and w2 are assumed to be shortened
int estimateBraidCancellation(size_t n, const Word& w1, const Word& w2, size_t sufficient_len = 100) {
  int result;
  if (w1.length() <= sufficient_len) {
    if (w2.length() <= sufficient_len) {
      return w1.length() + w2.length() - shortenBraid2(n, w1 * w2).length();
    } else {
      const auto s2 = w2.initialSegment(sufficient_len);
      return w1.length() + s2.length() - shortenBraid2(n, w1 * s2).length();
    }
  } else {
    const auto s1 = w1.terminalSegment(sufficient_len);
    if (w2.length() <= sufficient_len) {
      return s1.length() + w2.length() - shortenBraid2(n, s1 * w2).length();
    } else {
      const auto s2 = w2.initialSegment(sufficient_len);
      return s1.length() + s2.length() - shortenBraid2(n, s1 * s2).length();
    }
  }
}

template <typename Stabilizer>
Word getReplacement(int gen);

template <>
Word getReplacement<StabilizerSquare>(int gen) {
  return Word(-gen);
}

template <>
Word getReplacement<StabilizerDoubleSquare>(int gen) {
  return Word({-gen, -gen, -gen});
}

//! Returns a pair (i, w') where i is the position where the flip was made
//! and w' is obtained from w by replacing w[i] with w[i]^{-1} or with w[i]^{-3}.
template <typename Stabilizer>
std::vector<std::pair<size_t, Word>> availableFlips(size_t n, size_t a, size_t b, const Word& w) {
  static const crag::braidgroup::BraidHasher<FF> hasher(n);

  std::vector<std::pair<size_t, Word>> result;

  if (a > b) {
    std::swap(a, b);
  }

  // 1. Prepare initial list of words
  size_t pos = 0;
  Permutation perm1(n), perm2(n);
  std::set<braid_hash_t> flip_hash;

  for (const auto l : w) {
    const auto ind = std::abs(l);

    auto strand1 = perm2[ind - 1];
    auto strand2 = perm2[ind];

    perm1.change(strand1, strand2);
    perm2.change(ind - 1, ind);

    if (strand1 > strand2) {
      std::swap(strand1, strand2);
    }

    // available flip
    if (strand1 == a && strand2 == b) {
      auto w1 = w.subword(0, pos);
      w1 *= getReplacement<Stabilizer>(l);
      w1 *= w.subword(pos + 1, w.size());

      const auto h = hasher(w1);

      if (flip_hash.find(h) == flip_hash.end()) {
        flip_hash.insert(h);
        result.push_back(std::make_pair(pos, w1));
      }
    }

    pos++;
  }

  // 2. Apply shortenBraid2
  return parallel::map(
      result, [n](const std::pair<size_t, Word>& p) { return std::make_pair(p.first, shortenBraid2(n, p.second)); });
}

size_t totalLength(const std::vector<Word>& v) {
  size_t result = 0;

  for (const auto& w : v) {
    result += w.length();
  }

  return result;
}

map<braid_hash_t, pair<Word, size_t>>::iterator findBest(map<braid_hash_t, pair<Word, size_t>>& m) {
  size_t best_size = std::numeric_limits<size_t>::max();
  map<braid_hash_t, pair<Word, size_t>>::iterator result;
  for (auto it = m.begin(); it != m.end(); ++it) {
    const auto s = std::get<size_t>(it->second);
    if (best_size > s) {
      result = it;
      best_size = s;
    }
  }
  return result;
}

pair<vector<Word>, Word> easyDescend(size_t n, const vector<Word>& v) {
  vector<Word> conjugators;
  for (int i = 1; i <= n; ++i) {
    conjugators.push_back(Word(i));
    conjugators.push_back(Word(-i));
  }

  Word c;
  const auto sz = v.size();
  vector<Word> result = v;

  for (bool progress = true; progress;) {
    progress = false;
    for (const auto d : conjugators) {
      vector<Word> result2 = result;
      for (auto& w : result2) {
        w.push_front(-d);
        w.push_back(d);
      }

      if (totalLength(result2) < totalLength(result)) {
        result = result2;
        c.push_back(d);
        progress = true;
        break;
      }
    }
  }
  if (c.length() != 0) {
    cout << "ED: " << totalLength(v) - totalLength(result) << endl;
  }
  return make_pair(result, c);
}

vector<Word> multiplyVectorByWordsOnBothSides(size_t n, const std::vector<Word>& v, const Word& g1, const Word& g2) {
  vector<Word> result;
  for (const auto& w : v) {
    result.push_back(shortenBraid2(n, g1 * w * g2));
    // result.push_back(shortenBraid(n, g1 * w * g2));
  }
  return result;
}

vector<Word> multiplyVectorByWordsOnBothSides_noreduction(const std::vector<Word>& v, const Word& g1, const Word& g2) {
  vector<Word> result;
  for (const auto& w : v) {
    result.push_back(g1 * w * g2);
  }
  return result;
}

template <typename Stabilizer>
std::pair<Word, Word>
removeLeftCloakingElementsForPair(size_t n, size_t a, size_t b, const Word& sig1, const Word& sig2) {
  auto flips1 = availableFlips<Stabilizer>(n, a, b, sig1);
  auto flips2 = availableFlips<Stabilizer>(n, a, b, sig2);

  const auto cmp = [](const std::pair<size_t, Word>& lhs, const std::pair<size_t, Word>& rhs) {
    return lhs.first < rhs.first;
  };

  // consider only flips in the FIRST third of sig1
  const auto it1 = std::lower_bound(flips1.begin(), flips1.end(), std::make_pair(sig1.size() / 3, Word()), cmp);
  flips1.erase(it1, flips1.end());

  // consider only flips in the FIRST third of sig2
  const auto it2 = std::lower_bound(flips2.begin(), flips2.end(), std::make_pair(sig2.size() / 3, Word()), cmp);
  flips2.erase(it2, flips2.end());

  const auto flip_pairs_count = flips1.size() * flips2.size();

  assert(flip_pairs_count != 0);

  const auto unwrap = [&](size_t ij) { return std::make_pair(ij / flips2.size(), ij % flips2.size()); };

  const auto lengths = parallel::map<size_t>(flip_pairs_count, [&](size_t ij) {
    // unwrap indices
    size_t i, j;
    std::tie(i, j) = unwrap(ij);
    const auto& f1 = flips1[i];
    const auto& f2 = flips2[j];

    const auto prefix_length = std::max(f1.second.length() / 3, f2.second.length() / 3);

    const auto& s1 = flips1[i].second.subword(0, prefix_length);
    const auto& s2 = flips2[j].second.subword(0, prefix_length);

    return shortenBraid2(n, -s1 * s2).length();
  });

  const auto it = std::min_element(lengths.begin(), lengths.end());
  const auto idx = std::distance(lengths.begin(), it);

  // unwrap indices
  size_t i, j;
  std::tie(i, j) = unwrap(idx);

  return std::make_pair(flips1[i].second, flips2[j].second);
}

template <typename Stabilizer>
std::vector<Word> removeLeftCloakingElements(size_t n, size_t a, size_t b, const std::vector<Word>& sig) {
  assert(sig.size() % 2 == 0);

  const auto pairs_count = sig.size() / 2;

  std::vector<Word> result;
  result.reserve(sig.size());

  for (size_t i = 0; i < pairs_count; ++i) {
    const auto r = removeLeftCloakingElementsForPair<Stabilizer>(n, a, b, sig[2 * i], sig[2 * i + 1]);
    result.push_back(r.first);
    result.push_back(r.second);
  }

  return result;
}

template <typename Stabilizer>
std::pair<Word, Word>
removeRightCloakingElementsForPair(size_t n, size_t a, size_t b, const Word& sig1, const Word& sig2) {
  std::pair<Word, Word> result;
  auto flips1 = availableFlips<Stabilizer>(n, a, b, sig1);
  auto flips2 = availableFlips<Stabilizer>(n, a, b, sig2);

  auto cmp = [](const std::pair<size_t, Word>& lhs, const std::pair<size_t, Word>& rhs) {
    return lhs.first < rhs.first;
  };

  // consider only flips in the LAST third of sig1
  const auto it1 = std::upper_bound(flips1.begin(), flips1.end(), std::make_pair(sig1.size() / 3, Word()), cmp);
  flips1.erase(flips1.begin(), it1);

  // consider only flips in the LAST third of sig2
  const auto it2 = std::upper_bound(flips2.begin(), flips2.end(), std::make_pair(sig2.size() / 3, Word()), cmp);
  flips2.erase(flips2.begin(), it2);

  const auto flip_pairs_count = flips1.size() * flips2.size();

  assert(flip_pairs_count != 0);

  const auto unwrap = [&](size_t ij) { return std::make_pair(ij / flips2.size(), ij % flips2.size()); };

  const auto lengths = parallel::map<size_t>(flip_pairs_count, [&](size_t ij) {
    // unwrap indices
    size_t i, j;
    std::tie(i, j) = unwrap(ij);

    const auto s1_full_size = flips1[i].second.size();
    const auto s2_full_size = flips2[j].second.size();

    const auto suffix_length = std::max(s1_full_size / 3, s2_full_size / 3);

    const auto& s1 = flips1[i].second.subword(s1_full_size - suffix_length, s1_full_size);
    const auto& s2 = flips2[j].second.subword(s2_full_size - suffix_length, s2_full_size);

    return shortenBraid2(n, s1 * -s2).length();
  });

  const auto it = std::min_element(lengths.begin(), lengths.end());
  const auto idx = std::distance(lengths.begin(), it);

  // unwrap indices
  size_t i, j;
  std::tie(i, j) = unwrap(idx);

  return std::make_pair(flips1[i].second, flips2[j].second);
}

template <typename Stabilizer>
std::vector<Word> removeRightCloakingElements(size_t n, size_t a, size_t b, const vector<Word>& sig) {
  assert(sig.size() % 2 == 0);

  const auto pairs_count = sig.size() / 2;

  std::vector<Word> result;
  result.reserve(sig.size());

  for (size_t i = 0; i < pairs_count; ++i) {
    const auto r = removeRightCloakingElementsForPair<Stabilizer>(n, a, b, sig[2 * i], sig[2 * i + 1]);
    result.push_back(r.first);
    result.push_back(r.second);
  }

  return result;
}

//! Precise test for properly uncloaked signatures, only for diagnostics.
bool isUncloakedCorrectly(size_t n, const PrivateKey& private_key, const Signature& s, const Word& s_uncloaked) {
  return areEqualBraids(n, s_uncloaked, (-private_key.w1() * s.encodedMessageHash() * private_key.w2()));
}

//! Heuristic test for properly uncloaked signatures available for the attack.
bool areUncloakedCorrectly(
    size_t n, const Signature& s1, const Signature& s2, const Word& s1_uncloaked, const Word& s2_uncloaked) {
  static braidgroup::FastConjugacyChecker<FF> conjugacy_checker(n, 0);

  const auto lhs = s1_uncloaked * -s2_uncloaked;
  const auto rhs = s1.encodedMessageHash() * -s2.encodedMessageHash();

  const auto are_conjugate = !conjugacy_checker.areNotConjugate(lhs, rhs);

  return are_conjugate;
}

template <typename T, typename Obfuscator, typename Encoder, typename Stabilizer>
boost::optional<std::pair<Word, Word>> removeCloakingElements(
    const Protocol<T, Obfuscator, Encoder, Stabilizer>& p,
    const PrivateKey&,
    const PublicKey<T>& public_key,
    const Signature& s1,
    const Signature& s2) {
  const auto n = p.publicParameters().n();
  const auto a = p.publicParameters().a();
  const auto b = p.publicParameters().b();

  const auto sigma_w1 = public_key.w1Projection().permutation().inverse();

  const auto strand1 = sigma_w1[a];
  const auto strand2 = sigma_w1[b];

  const auto without_v1 =
      removeLeftCloakingElementsForPair<Stabilizer>(n, strand1, strand2, s1.signature(), s2.signature());

  //  std::cout << "Removed v1 from s1: " << std::boolalpha
  //            << areEqualBraids(n, without_v1.first, -s1.v1() * s1.signature()) << std::endl;
  //  std::cout << "Removed v1 from s2: " << std::boolalpha
  //            << areEqualBraids(n, without_v1.second, -s2.v1() * s2.signature()) << std::endl;

  const auto without_v1_v2 =
      removeRightCloakingElementsForPair<Stabilizer>(n, strand1, strand2, without_v1.first, without_v1.second);

  //  std::cout << "Removed v2 from s1: " << std::boolalpha
  //            << areEqualBraids(n, without_v1_v2.first, -s1.v1() * s1.signature() * -s1.v2()) << std::endl;
  //  std::cout << "Removed v2 from s2: " << std::boolalpha
  //            << areEqualBraids(n, without_v1_v2.second, -s2.v1() * s2.signature() * -s2.v2()) << std::endl;

  const auto cmp = [](const std::pair<size_t, Word>& lhs, const std::pair<size_t, Word>& rhs) {
    return lhs.second.size() < rhs.second.size();
  };

  const auto s1_flips = availableFlips<Stabilizer>(n, strand1, strand2, without_v1_v2.first);
  //  std::cout << "1st signature flips: " << s1_flips.size() << std::endl;
  const auto s1_uncloaked = std::min_element(s1_flips.begin(), s1_flips.end(), cmp)->second;

  //  std::cout << "Removed v from s1: " << std::boolalpha
  //            << areEqualBraids(n, s1_uncloaked, -private_key.w1() * s1.encodedMessageHash() * private_key.w2())
  //            << std::endl;

  const auto s2_flips = availableFlips<Stabilizer>(n, strand1, strand2, without_v1_v2.second);
  //  std::cout << "2nd signature flips: " << s2_flips.size() << std::endl;
  const auto s2_uncloaked = std::min_element(s2_flips.begin(), s2_flips.end(), cmp)->second;

  //  std::cout << "Removed v from s2: " << std::boolalpha
  //            << areEqualBraids(n, s2_uncloaked, -private_key.w1() * s2.encodedMessageHash() * private_key.w2())
  //            << std::endl;

  if (areUncloakedCorrectly(n, s1, s2, s1_uncloaked, s2_uncloaked)) {
    return std::make_pair(s1_uncloaked, s2_uncloaked);
  }

  return boost::none;
}

//! Returns minimal and maximal indices of generators in a vector of words.
std::pair<size_t, size_t> range(const vector<Word>& words) {
  size_t min_index = std::numeric_limits<size_t>::max();
  size_t max_index = 0;

  for (const auto& word : words) {
    for (const auto g : word) {
      const auto index_abs = static_cast<size_t>(std::abs(g));

      min_index = std::min(min_index, index_abs);
      max_index = std::max(max_index, index_abs);
    }
  }

  return std::make_pair(min_index, max_index);
}

//! Conjugates a word by small delta in degree = k.
Word conjugateBySmallDelta(const Word& word, size_t k) {
  std::vector<int> result;

  for (const auto g : word) {
    result.push_back(g > 0 ? g - k : g + k);
  }

  return Word(std::move(result));
}

//! Conjugates a tuple of words by small delta in degree = k.
vector<Word> conjugateBySmallDelta(const vector<Word>& words, size_t k) {
  vector<Word> result;
  result.reserve(words.size());

  for (const auto& w : words) {
    result.push_back(conjugateBySmallDelta(w, k));
  }

  return result;
}

vector<Word> conjugateByDelta(const vector<Word>& v, size_t max_index) {
  vector<Word> result;
  for (const auto& w : v) {
    Word new_w;
    for (const auto& l : w) {
      new_w.push_back(l > 0 ? max_index - l + 1 : -max_index - l - 1);
    }
    result.push_back(new_w);
  }
  return result;
}

boost::optional<vector<Word>> pushToLowerRank(size_t n, const vector<Word>& words) {
  const auto r = range(words);

  const auto left_handle_free_words = parallel::map(words, [&](const Word& w) {
    LinkedBraidStructure lbs(n - 1, w);
    lbs.removeLeftHandles();
    return lbs.translateIntoWord();
  });

  auto new_range = range(left_handle_free_words);

  if (new_range.first > r.first || new_range.second < r.second) {
    return left_handle_free_words;
  }

  const auto right_handle_free_words = parallel::map(words, [&](const Word& w) {
    LinkedBraidStructure lbs(n - 1, w);
    lbs.removeRightHandles();
    return lbs.translateIntoWord();
  });

  new_range = range(right_handle_free_words);

  if (new_range.first > r.first || new_range.second < r.second) {
    return right_handle_free_words;
  }

  return {};
}

void resetEnumeration(
    const vector<Word>& v,
    const Word& y,
    map<braid_hash_t, vector<Word>>& hash_values,
    map<braid_hash_t, pair<Word, size_t>>& checked_elts,
    map<braid_hash_t, pair<Word, size_t>>& unchecked_elts,
    const crag::braidgroup::BraidHasher<FF>& hasher) {
  hash_values.clear();
  checked_elts.clear();
  unchecked_elts.clear();
  const auto new_hash = hasher(v);
  hash_values[new_hash] = v;
  unchecked_elts[new_hash] = make_pair(y, totalLength(v));
}

// Returns fundamental braid for BKL form.
// Conjugation by it shifts indices: i -> i - 1.
Word smallDelta(size_t n) {
  std::vector<int> delta;
  delta.reserve(n);

  for (int i = 1; i < n; ++i) {
    delta.push_back(i);
  }

  return Word(std::move(delta));
}

std::vector<Word> constructConjugators(size_t n) {
  std::vector<Word> result;
  result.reserve(2 * n);


  for (int i = 1; i < n; ++i) {
    result.push_back(Word(i));
    result.push_back(Word(-i));
  }

  return result;
}

//! For a vector of words returns a vector of shortened conjugates.
std::vector<Word> conjugate(size_t n, const std::vector<Word>& words, const Word& c) {
  std::vector<Word> conjugated_words;
  conjugated_words.reserve(words.size());

  for (const auto& w : words) {
    conjugated_words.push_back(shortenBraid2(n, (-c) * w * c));
  }

  return conjugated_words;
}

//! For a vector of words returns a vector of conjugates.
std::vector<Word> conjugateWithoutShortening(const std::vector<Word>& words, const Word& c) {
  std::vector<Word> conjugated_words;
  conjugated_words.reserve(words.size());

  for (const auto& w : words) {
    conjugated_words.push_back((-c) * w * c);
  }

  return conjugated_words;
}

std::vector<std::vector<Word>> shortenTuplesParallel(size_t n, const std::vector<std::vector<Word>>& tuples) {
  if (tuples.empty()) {
    return tuples;
  }

  const auto tuples_count = tuples.size();
  const auto tuple_size = tuples.front().size();

  std::vector<std::vector<Word>> result(tuples_count, std::vector<Word>(tuple_size));

  const auto words_count = tuples_count * tuple_size;

  parallel::forEach(words_count, [&](size_t i) {
    const auto tuple_idx = i / tuple_size;
    const auto word_idx = i % tuple_size;

    result[tuple_idx][word_idx] = shortenBraid2(n, tuples[tuple_idx][word_idx]);
  });

  return result;
}

boost::optional<braid_hash_t> generateNewElts(
    size_t n,
    map<braid_hash_t, vector<Word>>& hash_values,
    map<braid_hash_t, vector<Word>>& hash_values2,
    map<braid_hash_t, pair<Word, size_t>>& checked_elts,
    map<braid_hash_t, pair<Word, size_t>>& unchecked_elts,
    const crag::braidgroup::BraidHasher<FF>& hasher,
    bool init_segments_as_conjugators = false) {
  // 1. Take the best unchecked instance and its characteristics
  const auto it = findBest(unchecked_elts);
  const auto cur_hash = it->first;
  const auto& cur_vec = hash_values[cur_hash];
  const auto cur_len = std::get<size_t>(it->second);
  const Word cur_y = std::get<Word>(it->second);
  checked_elts[cur_hash] = it->second;
  unchecked_elts.erase(it);

  const auto r = range(cur_vec);
  cout << "Current |y| = " << cur_len << ", [" << r.first << "," << r.second << "]" << endl;

  // 2. If cur_vec does not involve x1, then conjugate by delta and reduce all indices
  if (size_t delta = r.first - 1) {
    cout << "Reducing range..." << endl;
    const auto c = smallDelta(n).power(delta);
    resetEnumeration(
        conjugateBySmallDelta(cur_vec, r.first - 1), cur_y * c, hash_values, checked_elts, unchecked_elts, hasher);
    return boost::none;
  }

  if (r.second > 4) {
    if (const auto new_vec = pushToLowerRank(n, cur_vec)) {
      cout << "Pushing to lower rank..." << endl;
      resetEnumeration(new_vec.get(), cur_y, hash_values, checked_elts, unchecked_elts, hasher);
      return boost::none;
    }
  }

  // 3. Form a set of transformations
  vector<Word> generators;
  generators.reserve(2 * n + cur_vec.size());

  if (!init_segments_as_conjugators) {
    for (int i = 1; i < r.second + 1; ++i) {
      generators.push_back(Word(i));
      generators.push_back(Word(-i));
    }
  }

  const auto l = init_segments_as_conjugators ? 250 : 10;

  for (size_t i = 0; i < cur_vec.size(); ++i) {
    generators.push_back(cur_vec[i].subword(0, l));
  }

  // 4. Apply transformations
  //  const auto conjugated_words =
  //      parallel::map<Word, std::vector<Word>>(generators, [&](const Word& c) { return conjugate(n, cur_vec, c); });

  std::vector<std::vector<Word>> conjugated_words;
  conjugated_words.reserve(generators.size());

  for (const auto& c : generators) {
    conjugated_words.push_back(conjugateWithoutShortening(cur_vec, c));
  }

  conjugated_words = shortenTuplesParallel(n, conjugated_words);

  for (size_t i = 0; i < conjugated_words.size(); ++i) {
    const auto& conj = generators[i];
    const auto& new_vec = conjugated_words[i];
    const auto new_hash = hasher(new_vec);
    const auto new_vec2 = conjugateByDelta(new_vec, r.second);
    const auto new_hash2 = hasher(new_vec2);
    const auto new_cur_y = cur_y * conj;

    if (hash_values.find(new_hash) == hash_values.end() && hash_values.find(new_hash2) == hash_values.end()) {
      hash_values[new_hash] = new_vec;
      unchecked_elts[new_hash] = make_pair(new_cur_y, totalLength(new_vec));

      if (hash_values2.find(new_hash) != hash_values2.end()) {
        cout << "We've done it 1!!!" << endl;
        return new_hash;
      }

      if (hash_values2.find(new_hash2) != hash_values2.end()) {
        cout << "We've done it 2!!!" << endl;
        Word delta_word = Word(Permutation::getHalfTwistPermutation(r.second + 1).geodesicWord());
        hash_values[new_hash2] = new_vec2;
        unchecked_elts[new_hash2] = make_pair(new_cur_y * delta_word, totalLength(new_vec2));
        return new_hash2;
      }
    }
  }

  return {};
}

// Check if some components perform poorly (happens due to poorly removed cloaking elements)
bool drop_poor_performing_components(
    size_t n,
    const crag::braidgroup::BraidHasher<FF>& hasher,
    vector<Word>& vec1,
    vector<Word>& vec2,
    map<braid_hash_t, vector<Word>>& hash_values1,
    map<braid_hash_t, pair<Word, size_t>>& checked_elts1,
    map<braid_hash_t, pair<Word, size_t>>& unchecked_elts1,
    map<braid_hash_t, vector<Word>>& hash_values2,
    map<braid_hash_t, pair<Word, size_t>>& checked_elts2,
    map<braid_hash_t, pair<Word, size_t>>& unchecked_elts2) {
  const auto h1 = findBest(checked_elts1)->first;
  const auto h2 = findBest(checked_elts2)->first;
  auto v1 = hash_values1[h1];
  auto v2 = hash_values2[h2];
  auto t1 = checked_elts1[h1];
  auto t2 = checked_elts2[h2];
  vector<int> length_change;
  vector<int> length_diff;
  for (auto i = 0u; i < v1.size(); ++i) {
    const int d = vec1[i].length() - v1[i].length();
    length_change.push_back(d);
    int diff = (int) vec1[i].length() - (int) vec2[i].length();
    length_diff.push_back(abs(diff));
    cout << d << " | ";
  }
  cout << endl;
  const int min_elt = *std::min_element(length_change.begin(), length_change.end());
  const int max_elt = *std::max_element(length_diff.begin(), length_diff.end());
  for (int i = 0; i < v1.size(); ++i) {
    const auto d = vec1[i].length() - v1[i].length();
    int diff = abs((int) vec1[i].length() - (int) vec2[i].length());
    if (d == min_elt) {
      // if (diff == max_elt) {
      v1.erase(v1.begin() + i);
      v2.erase(v2.begin() + i);
      vec1.erase(vec1.begin() + i);
      vec2.erase(vec2.begin() + i);
      cout << "    Drop component #" << i << endl;
      resetEnumeration(v1, std::get<Word>(t1), hash_values1, checked_elts1, unchecked_elts1, hasher);
      resetEnumeration(v2, std::get<Word>(t2), hash_values2, checked_elts2, unchecked_elts2, hasher);
      return true;
    }
  }
  return false;
}

template <typename T, typename Obfuscator, typename Encoder, typename Stabilizer, typename URNG>
boost::optional<PrivateKey> attack(
    const Protocol<T, Obfuscator, Encoder, Stabilizer>& p,
    const PrivateKey& private_key,
    const PublicKey<T>& public_key,
    URNG& g) {
  const auto n = p.publicParameters().n();
  const auto hash_size = p.encoder().hashSize();

  static const crag::braidgroup::FastIdentityChecker<FF> checker(n);
  static const crag::braidgroup::BraidHasher<FF> hasher(n);

  const size_t signatures_pairs_count = 3;
  const size_t signatures_count = 2 * signatures_pairs_count;

  // S_i
  std::vector<Signature> signatures;
  signatures.reserve(signatures_count);

  // P_i
  std::vector<Word> reduced_signatures;
  reduced_signatures.reserve(signatures_count);

  size_t uncloaked_pairs = 0;
  size_t pair_index = 0;

  while (uncloaked_pairs < signatures_pairs_count) {
    std::cout << "Uncloaking signatures pair #" << pair_index++;

    const auto s1 = p.sign(walnut::randomMessageHash(hash_size, g), private_key, g);
    const auto s2 = p.sign(walnut::randomMessageHash(hash_size, g), private_key, g);

    std::cout << ", |sig| = " << s1.signature().size() << ", " << s2.signature().size() << ", result ";

    if (const auto uncloaked_sigs = removeCloakingElements(p, private_key, public_key, s1, s2)) {
      signatures.push_back(s1);
      signatures.push_back(s2);

      reduced_signatures.push_back(uncloaked_sigs->first);
      reduced_signatures.push_back(uncloaked_sigs->second);

      ++uncloaked_pairs;
      std::cout << "+" << std::endl;
    } else {
      std::cout << "X" << std::endl;
    }
  }

  // DIAGNOSTICS
  std::cout << "Cloaking elements success: ";
  for (const auto b : parallel::bmap(signatures_count, [&](size_t i) {
         return isUncloakedCorrectly(n, private_key, signatures[i], reduced_signatures[i]);
       })) {
    std::cout << (b ? "+" : "X");
  }
  std::cout << std::endl;

  // solve system of conjugacy equations
  std::vector<Word> vec_lhs;
  vec_lhs.reserve(signatures_count);

  std::vector<Word> vec_rhs;
  vec_rhs.reserve(signatures_count);

  for (size_t i = 1; i < signatures_count; ++i) {
    // P_{i-1} * P_i^{-1}
    vec_lhs.push_back(reduced_signatures[i - 1] * -reduced_signatures[i]);

    // E(H(m_{i-1})) * E(H(m_{i}))^{-1}
    vec_rhs.push_back(signatures[i - 1].encodedMessageHash() * -signatures[i].encodedMessageHash());
  }

  vec_lhs = parallel::map(vec_lhs, [&](const Word& w) { return shortenBraid2(n, w); });
  vec_rhs = parallel::map(vec_rhs, [&](const Word& w) { return shortenBraid2(n, w); });

  std::map<braid_hash_t, std::vector<Word>> hash_values1;
  std::map<braid_hash_t, std::pair<Word, size_t>> checked_elts1;
  std::map<braid_hash_t, std::pair<Word, size_t>> unchecked_elts1;
  const auto h1 = hasher(vec_lhs);
  hash_values1[h1] = vec_lhs;
  unchecked_elts1[h1] = std::make_pair(Word(), totalLength(vec_lhs));

  std::map<braid_hash_t, std::vector<Word>> hash_values2;
  std::map<braid_hash_t, std::pair<Word, size_t>> checked_elts2;
  std::map<braid_hash_t, std::pair<Word, size_t>> unchecked_elts2;
  const auto h2 = hasher(vec_rhs);
  hash_values2[h2] = vec_rhs;
  unchecked_elts2[h2] = std::make_pair(Word(), totalLength(vec_rhs));

  // One special iteration to drop the length
  generateNewElts(n, hash_values1, hash_values2, checked_elts1, unchecked_elts1, hasher, true);

  // 5. Start enumeration
  for (int step = 0; step < 200; ++step) {
    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);

    std::cout << ".................................................." << std::endl;
    std::cout << "Step #" << step << std::put_time(&tm, ",  tm = %H-%M-%S") << std::endl;
    std::pair<Word, size_t> t1, t2;

    if (const auto h = generateNewElts(n, hash_values1, hash_values2, checked_elts1, unchecked_elts1, hasher)) {
      const auto hash = h.get();
      t1 = unchecked_elts1[hash];
      t2 = unchecked_elts2.find(hash) != unchecked_elts2.end() ? unchecked_elts2[hash] : checked_elts2[hash];
    } else if (const auto h = generateNewElts(n, hash_values2, hash_values1, checked_elts2, unchecked_elts2, hasher)) {
      const auto hash = h.get();
      t1 = unchecked_elts1.find(hash) != unchecked_elts1.end() ? unchecked_elts1[hash] : checked_elts1[hash];
      t2 = unchecked_elts2[hash];
    } else {
      // Drop bad components. This eliminates signatures where we did poor job on removing cloaking elements
      if (step >= 5 && step % 5 == 0 && vec_lhs.size() > 3 && !checked_elts1.empty() && !checked_elts2.empty()) {
        drop_poor_performing_components(
            n,
            hasher,
            vec_lhs,
            vec_rhs,
            hash_values1,
            checked_elts1,
            unchecked_elts1,
            hash_values2,
            checked_elts2,
            unchecked_elts2);
      }
      continue;
    }

    // Check correctness of the obtained keys
    const Word w1_ = std::get<0>(t2) * -std::get<0>(t1);
    const Word w2_ = -signatures[0].encodedMessageHash() * w1_ * reduced_signatures[0];
    const auto h1 = hasher(multiplyVectorByWordsOnBothSides_noreduction(vec_rhs, -w1_, w1_));
    const auto h2 = hasher(vec_lhs);

    if (h1 != h2) {
      std::cout << "Conjugator is not correct!!!!" << std::endl;
      exit(1);
    } else {
      std::cout << "Conjugator seems to be correct" << std::endl;
    }

    return PrivateKey(w1_, w2_);
  }

  return boost::none;
}

template <typename T, typename Obfuscator, typename Encoder, typename Stabilizer, typename URNG>
bool checkFakePrivateKey(
    const Protocol<T, Obfuscator, Encoder, Stabilizer>& p,
    const PrivateKey& fake_private_key,
    const PrivateKey& original_private_key,
    const PublicKey<T>& original_public_key,
    URNG& g) {
  const auto hash_size = p.encoder().hashSize();
  const auto signatures_count = 10;

  std::vector<Signature> fake_signatures;
  fake_signatures.reserve(signatures_count);

  std::vector<Signature> original_signatures;
  original_signatures.reserve(signatures_count);

  for (size_t i = 0; i < signatures_count; ++i) {
    const auto random_msg_hash = randomMessageHash(hash_size, g);
    const auto fake_signature = p.signWithoutCloakingEls(random_msg_hash, fake_private_key);
    const auto original_signature = p.signWithoutCloakingEls(random_msg_hash, original_private_key);

    fake_signatures.push_back(fake_signature);
    original_signatures.push_back(original_signature);
  }

  const auto n = p.publicParameters().n();

  const auto short_fake_signatures = parallel::map<Signature, Word>(
      fake_signatures, [n](const Signature& s) { return shortenBraid(n, s.signature()); });

  const auto short_original_signatures = parallel::map<Signature, Word>(
      original_signatures, [n](const Signature& s) { return shortenBraid(n, s.signature()); });

  std::ostringstream sig_cmp;
  sig_cmp << "Signature comparison fake vs real: ";

  std::ostringstream sig_length;
  sig_length << "(|fake sig|, |real sig|) = ";

  bool same_signatures = true;

  for (size_t i = 0; i < signatures_count; ++i) {
    const auto fake_signature = short_fake_signatures[i];
    const auto original_signature = short_original_signatures[i];

    LinkedBraidStructure lbs(n - 1, -fake_signature * original_signature);
    lbs.removeLeftHandles();

    same_signatures &= (lbs.size() == 0);

    sig_cmp << ((lbs.size() == 0) ? "+" : "X");
    sig_length << "(" << fake_signature.size() << ", " << original_signature.size() << "), ";
  }

  sig_cmp << std::endl;
  sig_length << std::endl;

  std::cout << sig_length.str() << sig_cmp.str() << "Same signatures = " << std::boolalpha << same_signatures
            << std::endl;

  // signatures are literally the same, so don't need to validate them
  if (same_signatures) {
    return true;
  }

  return std::all_of(fake_signatures.begin(), fake_signatures.end(), [&](const Signature& s) {
    return p.verify(s, original_public_key);
  });
}

using GF32 = finitefield::FieldElement<finitefield::IdealGeneratedByPolynomial<finitefield::ZZ<2>, 1, 0, 1, 0, 0, 1>>;
using GF256 =
    finitefield::FieldElement<finitefield::IdealGeneratedByPolynomial<finitefield::ZZ<2>, 1, 1, 0, 1, 1, 0, 0, 0, 1>>;

template <typename Stabilizer = StabilizerSquare>
auto getProtocolFor128BitsSecurity(size_t seed) {
  const size_t n = 8;
  const size_t hash_size = 256;
  const size_t L = 15;

  const auto random_parameters = randomParameters<GF32, Stabilizer>(n, seed).wMinLength(132).wMaxLength(150);
  const auto random_encoder = randomEncoder(n, hash_size, seed);
  //  const auto random_encoder = randomAdvancedEncoder(n, n, hash_size, seed);

  // corresponding to parameters generator of elements from a stabilizer
  const auto stabilizer = Stabilizer(random_parameters, L);

  return walnut::getProtocol(random_parameters, random_encoder, stabilizer);
}

template <typename Stabilizer = StabilizerSquare>
auto getProtocolFor256BitsSecurity(size_t seed) {
  const size_t n = 8;
  const size_t hash_size = 512;
  const size_t L = 30;

  const auto random_parameters = randomParameters<GF256, Stabilizer>(n, seed).wMinLength(287).wMaxLength(300);
  const auto random_encoder = randomEncoder(n, hash_size, seed);
  //  const auto random_encoder = randomAdvancedEncoder(n, n, hash_size, seed);

  // corresponding to parameters generator of elements from a stabilizer
  const auto stabilizer = Stabilizer(random_parameters, L);

  return walnut::getProtocol(random_parameters, random_encoder, stabilizer);
}

template <typename Stabilizer = StabilizerDoubleSquare>
auto getProtocolFor256BitsSecurityN11(size_t seed) {
  const size_t n = 11;
  const size_t hash_size = 512;
  const size_t L = 30;

  const auto random_parameters = randomParameters<GF256, Stabilizer>(n, seed).wMinLength(287).wMaxLength(300);
  const auto random_encoder = randomEncoder(n, hash_size, seed);
  //  const auto random_encoder = randomAdvancedEncoder(n, n, hash_size, seed);

  // corresponding to parameters generator of elements from a stabilizer
  const auto stabilizer = Stabilizer(random_parameters, L);

  return walnut::getProtocol(random_parameters, random_encoder, stabilizer);
}
} // namespace walnut
} // namespace crag

#endif // CRAG_WALNUT_ATTACK_H
