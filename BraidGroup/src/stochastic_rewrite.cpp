#include "stochastic_rewrite.h"
#include "braid_group.h"
#include <numeric>

namespace crag {
namespace stochasticrewrite {

std::vector<size_t> calculateRs(const std::vector<size_t>& p, size_t init) {
  std::vector<size_t> result;

  result.push_back(init);

  for (size_t i = 0; i < p.size(); ++i) {
    result.push_back(result.back() + p[i]);
  }

  return result;
}


std::vector<Word> calculateYGensInBGens(const std::vector<size_t>& p) {
  const auto rs = calculateRs(p);

  std::vector<Word> result;

  for (size_t i = 1; i < rs.size(); ++i) {
    for (size_t j = rs[i - 1]; j < rs[i]; ++j) {
      std::vector<int> w(rs[i] - j);
      std::iota(w.begin(), w.end(), j);
      result.push_back(Word(w));
    }
  }

  return result;
}


std::vector<Word> calculateBGensInYGens(const std::vector<size_t>& p) {
  const auto rs = calculateRs(p);

  std::vector<Word> result;

  for (size_t i = 1; i < rs.size(); ++i) {
    for (int j = static_cast<int>(rs[i - 1]); j < static_cast<int>(rs[i]) - 1; ++j) {
      result.push_back(Word({j, -(j + 1)}));
    }

    result.push_back(Word({static_cast<int>(rs[i]) - 1}));
  }

  return result;
}


std::vector<Word> calculateYRelations(const std::vector<size_t>& p) {
  std::vector<Word> result;

  const auto rs = calculateRs(p);
  const auto n = rs.back();

  const auto rels = braidgroup::getRelations(n);

  const auto b_gens = calculateBGensInYGens(p);

  for (const auto& r : rels) {
    result.push_back(r.replaceGenerators(b_gens));
  }

  return result;
}


std::vector<Word> calculateAdditionalRelations(const std::vector<size_t>& p) {
  const auto rs = calculateRs(p);

  std::vector<Word> result;

  for (size_t k = 1; k < rs.size(); ++k) {
    for (int i = rs[k - 1]; i < static_cast<int>(rs[k]); ++i) {
      for (int j = i + 1; j <= static_cast<int>(rs[k]) - 1; ++j) {
        result.push_back(Word({j, i, static_cast<int>(rs[k]) - 1, 1 - j, -i}));
      }
    }
  }

  return result;
}


std::vector<Word> calculateSRelations(const std::vector<size_t>& p) {
  std::vector<Word> result;

  const auto y_relations = calculateYRelations(p);
  const auto additional_relations = calculateAdditionalRelations(p);

  auto relations = y_relations;
  relations.insert(relations.end(), additional_relations.begin(), additional_relations.end());

  for (const auto& r : relations) {
    result.push_back(r);
    for (size_t i = 1; i < r.length(); ++i) {
      auto t = result.back();
      t.cyclicLeftShift();
      result.push_back(t);
    }
  }

  return result;
}


std::multimap<Word, Word> calculateRewritingRules(const std::vector<size_t>& p) {
  std::multimap<Word, Word> result;

  const auto rels = calculateSRelations(p);

  for (const auto& r : rels) {
    result.emplace(r.segment(0, 2), r.terminalSegment(2) ^ -1);
    result.emplace(r.segment(0, 2) ^ -1, r.terminalSegment(2));
  }

  return result;
}

} // namespace stochasticrewrite
} // namespace crag
