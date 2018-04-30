#include "GCD.h"

namespace crag {
namespace finitefield {

template <>
std::pair<int, int> getNormed<int>(const int& x) {
  if (x < 0) {
    return std::make_pair(-x, -1);
  }

  return std::make_pair(x, 1);
}
} // namespace finitefield
} // namespace crag