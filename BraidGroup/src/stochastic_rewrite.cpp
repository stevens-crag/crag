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


std::vector<std::vector<int>> calculateYGensInBGens(const std::vector<size_t>& p) {
  const auto rs = calculateRs(p);

  std::vector<std::vector<int>> result;

  for (size_t i = 1; i < rs.size(); ++i) {
    for (size_t j = rs[i - 1]; j < rs[i]; ++j) {
      std::vector<int> w(rs[i] - j);
      std::iota(w.begin(), w.end(), j);
      result.push_back(w);
    }
  }

  return result;
}


std::vector<std::vector<int>> calculateBGensInYGens(const std::vector<size_t>& p) {
  const auto rs = calculateRs(p);

  std::vector<std::vector<int>> result;

  for (size_t i = 1; i < rs.size(); ++i) {
    for (int j = static_cast<int>(rs[i - 1]); j < static_cast<int>(rs[i]) - 1; ++j) {
      result.push_back({j, -(j + 1)});
    }

    result.push_back({static_cast<int>(rs[i]) - 1});
  }

  return result;
}


std::vector<std::vector<int>> calculateYRelations(const std::vector<size_t>& p) {
  std::vector<std::vector<int>> result;

  const auto rs = calculateRs(p);
  const auto n = rs.back();

  const auto rels = braidgroup::getRelations(n);

  const auto b_gens = calculateBGensInYGens(p);

  for (const auto& r : rels) {
    result.push_back(replaceGenerators(r.toVector(), b_gens));
  }

  return result;
}


std::vector<std::vector<int>> calculateAdditionalRelations(const std::vector<size_t>& p) {
  const auto rs = calculateRs(p);

  std::vector<std::vector<int>> result;

  for (size_t k = 1; k < rs.size(); ++k) {
    for (int i = rs[k - 1]; i < static_cast<int>(rs[k]); ++i) {
      for (int j = i + 1; j <= static_cast<int>(rs[k]) - 1; ++j) {
        result.push_back({j, i, static_cast<int>(rs[k]) - 1, 1 - j, -i});
      }
    }
  }

  return result;
}


std::vector<std::vector<int>> calculateSRelations(const std::vector<size_t>& p) {
  std::vector<std::vector<int>> result;

  const auto y_relations = calculateYRelations(p);
  const auto additional_relations = calculateAdditionalRelations(p);

  auto relations = y_relations;
  relations.insert(relations.end(), additional_relations.begin(), additional_relations.end());

  for (const auto& r : relations) {
    result.push_back(r);
    for (size_t i = 1; i < r.size(); ++i) {
      auto t = result.back();
      result.push_back(cyclicLeftShift(result.back()));
    }
  }

  return result;
}


std::multimap<std::vector<int>, std::vector<int>> calculateRewritingRules(const std::vector<size_t>& p) {
  std::multimap<std::vector<int>, std::vector<int>> result;

  const auto rels = calculateSRelations(p);

  for (const auto& r : rels) {
    std::vector<int> beg{r[0], r[1]};
    std::vector<int> rest(r.begin() + 2, r.end());

    result.emplace(beg, invert(rest));
    result.emplace(invert(beg), rest);
  }

  return result;
}


std::vector<int> replaceGenerators(const std::vector<int>& w, const std::vector<std::vector<int>>& images) {
  std::vector<int> result;

  std::vector<std::vector<int>> inv_images;

  for (const auto& image : images) {
    inv_images.push_back(invert(image));
  }

  for (auto a : w) {
    if (a > 0) {
      append(result, images[a - 1].begin(), images[a - 1].end());
    } else {
      append(result, inv_images[-a - 1].begin(), inv_images[-a - 1].end());
    }
  }

  return result;
}


std::vector<int> invert(const std::vector<int>& w) {
  std::vector<int> result;

  for (int i = w.size() - 1; i >= 0; --i) {
    result.push_back(-w[i]);
  }

  return result;
}


std::vector<int> cyclicLeftShift(const std::vector<int>& w) {
  auto result = w;
  
  std::rotate(result.begin(), result.begin() + 1, result.end());

  return result;
}

} // namespace stochasticrewrite
} // namespace crag
