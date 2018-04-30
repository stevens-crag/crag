#include "parallel.h"

#include <algorithm>

namespace crag {
namespace parallel {

size_t getHardwareConcurrency() {
  return std::max(1u, std::thread::hardware_concurrency());
}
}
}
