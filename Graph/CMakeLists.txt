#Specify the version being used aswell as the language
cmake_minimum_required(VERSION 2.8.12)
cmake_policy(VERSION 2.8.12)

include("../cmake/common.cmake")

crag_library(Graph
  FSA
  FSARep
  Graph
  GraphConcept
  GraphDrawingAttributes
  GraphRep
  GraphType
  RandomFSA
)

target_link_libraries(Graph
  PRIVATE Alphabet
  PUBLIC crag_general
  PRIVATE gmp
  PUBLIC ranlib
)

crag_main(test_rand Graph)
