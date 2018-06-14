#include "kayawood.h"

namespace crag {
namespace kayawood {

StochasticRewriteObfuscator getStochasticRewriteObfuscator(size_t n, size_t seed) {
  std::mt19937_64 g(seed);
  return StochasticRewriteObfuscator(stochasticrewrite::randomPartition(n - 1, 3, n - 1, g), 5, 10, 3, seed);
}
} // namespace kayawood
} // namespace crag
