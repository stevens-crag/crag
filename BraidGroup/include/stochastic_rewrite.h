#pragma once

#ifndef CRAG_STOCHASTIC_REWRITE_H
#define CRAG_STOCHASTIC_REWRITE_H

#include <boost/random/uniform_int_distribution.hpp>
#include <map>
#include <vector>

#include "Word.h"

namespace crag {
namespace stochasticrewrite {

/// Returns p1, ..., pl s. t. p1 + ... pl = n, pi >= min_value, pi <= max_value.
/// TODO: now the distribution is not uniform.
template <typename URNG>
std::vector<size_t> randomPartition(size_t n, size_t min_value, size_t max_value, URNG& g) {
  static boost::random::uniform_int_distribution<int> dist01(0, 1);

  std::vector<size_t> result;

  size_t p = 0;
  for (auto i = 0u; i < n - min_value; ++i) {
    ++p;

    if ((p >= min_value && dist01(g) == 0) || p == max_value) {
      result.push_back(p);
      p = 0;
    }
  }

  result.push_back(p + min_value);

  return result;
}

/// Returns r1 = init, r2 = r1 + p1, ..., ri = r[i-1] + pi, ..., rl = n.
std::vector<size_t> calculateRs(const std::vector<size_t>& p, size_t init = 1);

/// Returns y-generators as words in b-gens (see the paper). p is a partition, n is the number of strands.
std::vector<Word> calculateYGensInBGens(const std::vector<size_t>& p);

/// Returns b-gens as words in y-gens (see the paper). p is a partition, n is the number of strands.
std::vector<Word> calculateBGensInYGens(const std::vector<size_t>& p);

/// Expreses w as a word in y_gens;
Word expressInYGens(const Word& w, const std::vector<Word>& b_gens);

/// Returns the set yRel(P).
std::vector<Word> calculateYRelations(const std::vector<size_t>& p);

/// Returns the additional relations.
std::vector<Word> calculateAdditionalRelations(const std::vector<size_t>& p);

/// Returns SRel(p).
std::vector<Word> calculateSRelations(const std::vector<size_t>& p);


/// Returns rewriting rules u -> {w1, ..., wn}. where |u| == 2.
std::multimap<Word, Word> calculateRewritingRules(const std::vector<size_t>& p);


/// Returns beginnings of subwords us, |u| == 2.
template <typename URNG>
std::vector<size_t> selectSubwords(const std::vector<size_t>& blocks, URNG& g) {
  const auto rs = calculateRs(blocks, 0);

  std::vector<size_t> result;

  for (size_t i = 0; i < rs.size() - 1; ++i) {
    boost::random::uniform_int_distribution<size_t> dist(rs[i], rs[i + 1] - 2);
    result.push_back(dist(g));
  }

  return result;
}


// Performs one iteration of rewriting process.
template <typename URNG>
Word rewrite(
    const Word& w, size_t min_block_size, size_t max_block_size, const std::multimap<Word, Word>& rules, URNG& g) {
  const auto blocks = randomPartition(w.length(), min_block_size, max_block_size, g);
  const auto subwords = selectSubwords(blocks, g);

  std::vector<int> v = w.toVector();
  std::vector<int> result;

  size_t segment_begin = 0;

  for (size_t i = 0; i < subwords.size(); ++i) {
    const auto subword = Word({v[subwords[i]], v[subwords[i] + 1]});
    const auto appropriate_rules_count = rules.count(subword);
    if (appropriate_rules_count > 0) {
      boost::random::uniform_int_distribution<int> dist(0, appropriate_rules_count - 1);

      auto it = rules.lower_bound(subword);
      std::advance(it, dist(g));
      const auto new_subword = it->second;

      result.insert(result.end(), v.begin() + segment_begin, v.begin() + subwords[i]);
      result.insert(result.end(), new_subword.begin(), new_subword.end());

      segment_begin = subwords[i] + 2;
    }
  }

  result.insert(result.end(), v.begin() + segment_begin, v.end());

  return Word(result).freelyReduce();
}


// Stochastic rewriting.
template <typename URNG>
Word stochasticRewrite(
    const Word& w, const std::vector<size_t>& partition, size_t min_block_size, size_t max_block_size, size_t iter_num,
    URNG& g) {
  const auto b_gens = calculateBGensInYGens(partition);

  auto w_y = w.replaceGenerators(b_gens);

  const auto rules = calculateRewritingRules(partition);

  for (size_t i = 0; i < iter_num; ++i) {
    w_y = rewrite(w_y, min_block_size, max_block_size, rules, g);
  }

  const auto y_gens = calculateYGensInBGens(partition);

  return w_y.replaceGenerators(y_gens);
}

} // namespace stochasticrewrite
} // namespace crag

#endif // CRAG_STOCHASTIC_REWRITE_H
