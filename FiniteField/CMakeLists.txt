#Specify the version being used as well as the language
cmake_minimum_required(VERSION 2.8.12)
cmake_policy(VERSION 2.8.12)

include("../cmake/common.cmake")

crag_library(FiniteField
  FiniteField
)

target_link_libraries(FiniteField
  PUBLIC Boost::boost
)

crag_test(TestFiniteField FiniteField)

