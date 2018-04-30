cmake_minimum_required(VERSION 3.8)

include("../cmake/common.cmake")

crag_library(Walnut
  walnut
  walnut_encoding
)

target_link_libraries(Walnut
  PUBLIC BraidGroup
  PUBLIC crag_general
  PUBLIC Elt
  PUBLIC Random
)

crag_main(walnut_main Walnut)

crag_test(test_walnut Walnut)
crag_test(test_walnut_attack Walnut)
