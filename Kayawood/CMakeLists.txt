cmake_minimum_required(VERSION 3.8)

include("../cmake/common.cmake")

crag_library(Kayawood
  kayawood
)

target_link_libraries(Kayawood
  PUBLIC general
  PUBLIC BraidGroup
  PUBLIC Elt
  PUBLIC Random
)

crag_main(kayawood_main Kayawood)

crag_test(test_kayawood Kayawood)
crag_test(test_kayawood_attack Kayawood)
