#Specify the version being used aswell as the language
cmake_minimum_required(VERSION 2.8.12)
cmake_policy(VERSION 2.8.12)

include("../cmake/common.cmake")

crag_library(Elt
  Word
  WordIterator
  WordRep
)

target_link_libraries(Elt
  PUBLIC Alphabet
  PUBLIC crag_general
  PRIVATE ranlib
)

crag_main(test_word Elt)
crag_main(ac_enum Elt)
