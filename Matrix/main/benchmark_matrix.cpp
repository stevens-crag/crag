#include <random>

#include <benchmark/benchmark.h>

#include "matrix.h"

static void BM_MatrixMultiplication(benchmark::State& state) {
  std::mt19937 g(1233);
  std::uniform_int_distribution<> d(-10, 10);

  size_t n = state.range(0);
  
  while (state.KeepRunning()) { 
    crag::matrix::Matrix<int> a(state.range(0));
    crag::matrix::Matrix<int> b(state.range(0));

    for (size_t i = 0; i < n; ++i) {
      for (size_t j = 0; j < n; ++j) {
        a(i, j) = d(g);
        b(i, j) = d(g);
      }
    }

    benchmark::DoNotOptimize(a * b);
  }

  state.SetComplexityN(state.range(0));
}


BENCHMARK(BM_MatrixMultiplication)->RangeMultiplier(2)->Range(1, 512)->Complexity();

BENCHMARK_MAIN();
