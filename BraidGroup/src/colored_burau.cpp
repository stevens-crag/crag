#include "colored_burau.h"

namespace crag {
namespace coloredburau {

Permutation permutation(size_t n, const Word& w) {
  Permutation p(n);

  for (auto it = w.rbegin(); it != w.rend(); ++it) {
    const auto abs_i = static_cast<size_t>(std::abs(*it));

    if ((abs_i == 0) || (abs_i >= n)) {
      throw std::invalid_argument("Index of a generator is out of range.");
    }

    p.change(abs_i - 1, abs_i);
  }

  return p;
}
}
}
