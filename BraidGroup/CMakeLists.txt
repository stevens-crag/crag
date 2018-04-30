cmake_minimum_required(VERSION 3.8)

include("../cmake/common.cmake")

crag_library(BraidGroup
  braid_group
  ThRightNormalForm
  ThRightNormalForm_uss
  ThRightNormalForm_sss
  LinkedBraidStructure
  DehornoyForm
  ShortBraidForm
  ThLeftNormalForm
  colored_burau
  stochastic_rewrite
)

target_link_libraries(BraidGroup
  PUBLIC Boost::boost
  PUBLIC crag_general
  PUBLIC Elt
  PUBLIC FiniteField
  PUBLIC Polynomial
  PUBLIC Matrix
  PUBLIC Random
)

crag_main(test_rightNF BraidGroup)
crag_main(test_leftNF BraidGroup)
crag_main(test_deh_form BraidGroup)

# crag_main(mainParser Alphabet Elt)

crag_test(test_cloaking_element BraidGroup)
crag_test(test_colored_burau BraidGroup)
crag_test(test_fast_identity_check BraidGroup)
crag_test(test_stochastic_rewrite BraidGroup)
crag_test(test_linked_braid_structure BraidGroup)
crag_test(test_fast_conjugacy_check BraidGroup)
