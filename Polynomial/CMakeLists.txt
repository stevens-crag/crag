cmake_minimum_required(VERSION 3.8)

include("../cmake/common.cmake")

crag_library(Polynomial
  polynomial
)

target_link_libraries(Polynomial
  PUBLIC Boost::boost
  PUBLIC FiniteField
)

crag_test(test_polynomial Polynomial)
