#include <chrono>

#include "slp.h"
#include "EndomorphismSLP.h"
#include "gmp_boost_pool_allocator.h"

typedef crag::slp::TVertexHashAlgorithms<
    crag::slp::hashers::SinglePowerHash,
    crag::slp::hashers::PermutationHash<crag::Permutation16>
> VertexHashAlgorithms;


int main() {
  gmp_pool_setup();
  CONSTEXPR_OR_CONST size_t RANK = 6;
  CONSTEXPR_OR_CONST size_t ENDOMORPHISMS_NUMBER = 2200;
  size_t seed = 112233;
  crag::UniformAutomorphismSLPGenerator<> generator(RANK, seed);

  typedef crag::slp::TVertexHashAlgorithms<
      crag::slp::hashers::SinglePowerHash,
      crag::slp::hashers::PermutationHash<crag::Permutation16>
  > VertexHashAlgorithms;

  auto slp = crag::EndomorphismSLP::composition(ENDOMORPHISMS_NUMBER, generator).image(1);

  VertexHashAlgorithms::Cache calculated_hashes;
  std::unordered_map<crag::slp::Vertex, crag::slp::Vertex> reduced_vertices;

  auto begin = std::chrono::system_clock::now();
  auto reduced = VertexHashAlgorithms::reduce_narrow_slp(slp, &calculated_hashes, &reduced_vertices);
  auto end = std::chrono::system_clock::now();

  std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << std::endl;

  return 0;
}


