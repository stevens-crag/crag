#include "braid_group.h"

namespace crag {
namespace braidgroup {

Word BraidGroup::twist(const Word& w) const {
  std::list<int> result;

  for (auto w_it = w.begin(); w_it != w.end(); ++w_it) {
    int g = *w_it;

    result.push_back(g < 0 ? (-g - n_) : (n_ - g));
  }

  return Word(std::move(result));
}

std::vector<Word> getRelations(size_t n) {
  if (n < 2) {
    throw std::invalid_argument("Expect at least 2 generators for braid group.");
  }

  std::vector<Word> relations;
  relations.reserve(((n - 1) * (n - 2)) / 2);

  for (int i = 1; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      if (j - i == 1) {
        relations.push_back(Word({i, j, i, -j, -i, -j}));
      } else {
        relations.push_back(Word({i, j, -i, -j}));
      }
    }
  }

  return relations;
}

Word getPureBraidSubgroupGenerator(size_t n, size_t i, size_t j) {
  if (!((1 <= i) && (i < j) && (j <= n))) {
    throw std::invalid_argument("Invalid indices for a pure braid subgroup generator.");
  }

  Word w;

  for (int k = 0; k < j - i - 1; ++k) {
    w.push_back(j - k - 1);
  }

  return w * Word(i) * Word(i) * -w;
}
} // namespace braidgroup
} // namespace crag
