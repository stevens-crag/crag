/*
 * profile_matching_new.cpp
 *
 *  Created on: Apr 2, 2013
 *      Author: dpantele
 */

#include <memory>
#include <functional>
#include <utility>
#include <vector>
#include <tuple>
#include <algorithm>
#include <chrono>

#include "slp.h"
#include "EndomorphismSLP.h"
//#include "gmp_boost_pool_allocator.h"

using crag::slp::Vertex;
using crag::slp::PreorderInspector;
typedef crag::slp::TerminalVertexTemplate<int> TerminalVertex;
using crag::slp::NonterminalVertex;
using crag::slp::MatchingTable;
using crag::UniformAutomorphismSLPGenerator;
using crag::EndomorphismSLP;

int main() {
//  gmp_pool_setup();
  int REPEAT = 50;
  const size_t RANK = 3;
  const size_t ENDOMORPHISMS_NUMBER = 100;
  size_t seed = 112233;
  UniformAutomorphismSLPGenerator<int> generator(RANK, seed);
  auto begin = std::chrono::system_clock::now();
  int count = REPEAT;
  while (--count >= 0) {
    auto image = EndomorphismSLP<int>::composition(ENDOMORPHISMS_NUMBER, generator).slp(1);

    Vertex reduced = reduce(image);
    auto end = std::chrono::system_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << std::endl;
  }

  auto end = std::chrono::system_clock::now();
  std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()/REPEAT << std::endl;
}


