#Specify the version being used aswell as the language
cmake_minimum_required(VERSION 2.8.12)
cmake_policy(VERSION 2.8.12)

include("../cmake/common.cmake")

crag_library(CryptoShiftConj
  ShftConjKeyGeneraion
  ShftConjKeyGenerationGarside
)

target_link_libraries(CryptoShiftConj
  PUBLIC BraidGroup
  PRIVATE ranlib
)

crag_main(attack_shifted_CP_uss CryptoShiftConj ranlib)
