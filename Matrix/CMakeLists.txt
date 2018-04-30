cmake_minimum_required(VERSION 3.8)

include("../cmake/common.cmake")

crag_library(Matrix
  matrix
)

crag_test(test_matrix Matrix)

crag_main(benchmark_matrix Matrix benchmark::benchmark)
