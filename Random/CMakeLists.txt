cmake_minimum_required(VERSION 3.8)

include("../cmake/common.cmake")

crag_library(Random
  random_word
)

target_link_libraries(Random
  PUBLIC Boost::boost
  PUBLIC crag_general
  PUBLIC Elt
)

crag_test(test_random_permutation Random)
crag_test(test_random_subset Random)
crag_test(test_random_word Random)
crag_test(test_shuffle Random)
