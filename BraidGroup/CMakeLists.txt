#Specify the version being used aswell as the language
cmake_minimum_required(VERSION 2.8.12)
cmake_policy(VERSION 2.8.12)

include("../cmake/common.cmake")

crag_library(BraidGroup
  BraidGroup
  ThRightNormalForm
  ThRightNormalForm_uss
  ThRightNormalForm_sss
  LinkedBraidStructure
  DehornoyForm
  ShortBraidForm
  ThLeftNormalForm
)

target_link_libraries(BraidGroup
  PUBLIC Elt
)

crag_main(test_rightNF BraidGroup)
crag_main(test_leftNF BraidGroup)
crag_main(test_deh_form BraidGroup)

# crag_main(mainParser Alphabet Elt)
