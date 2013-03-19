
/*
 * endomorphism_slp_tests.cpp
 *
 *  Created on: Mar 19, 2013
 *      Author: pmorar
 */
#include <iostream>
#include <chrono>
#include "EndomorphismSLP.h"

typedef std::chrono::high_resolution_clock our_clock;
using namespace crag;

void composition_speed_test() {
  std::cout << "Time to compose (num iterations, num of composed elements, rank)" << std::endl;
  for (auto rank : {3, 5, 10, 20}) {
    UniformAutomorphismSLPGenerator<int> rnd(rank);
    for (auto size : {1000, 2000}) {
      auto start_time = our_clock::now();
      const int num = 10;
      for (int i = 0; i < num; ++i)
        EndomorphismSLP<int>::composition(size, rnd);
      auto time = our_clock::now() - start_time;
      auto time_in_ms = std::chrono::duration_cast<std::chrono::milliseconds>(time);
      std::cout << "(num=" << num << ", size=" << size << ", rank=" << rank << "): " <<
                   time_in_ms.count() << "ms, " << time_in_ms.count() / num << "ms per iteration" << std::endl;
    }
  }
}

void composition_statistics_test() {
  int heights_sum = 0;
  int heights_max = 0;
  for (auto rank : {3, 5, 10, 20}) {
    crag::UniformAutomorphismSLPGenerator<int> rnd(rank);
    for (auto size : {1000, 2000}) {
      auto start_time = our_clock::now();
      const int num = 10;
      for (int i = 0; i < num; ++i)
        crag::EndomorphismSLP<int>::composition(size, rnd);
      auto time = our_clock::now() - start_time;
      auto time_in_ms = std::chrono::duration_cast<std::chrono::milliseconds>(time);
      std::cout << "(num=" << num << ", size=" << size << ", rank=" << rank << "): " <<
                   time_in_ms.count() << "ms, " << time_in_ms.count() / num << "ms per iteration" << std::endl;
    }
  }
}

int main() {
  composition_speed_test();
}
