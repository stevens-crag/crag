#Specify the version being used aswell as the language
cmake_minimum_required(VERSION 2.8.12)
cmake_policy(VERSION 2.8.12)

include("../cmake/common.cmake")

crag_library(CryptoAGG
  AAGKeyGeneration
  AAGChallengeGeneration
  LengthAttack_A1
  LengthAttack_A2
  LengthAttack_A3
  MajorDump
)

target_link_libraries(CryptoAGG
  PUBLIC BraidGroup
  PRIVATE ranlib
)

crag_main(test_length_attack CryptoAGG ranlib)
crag_main(test_sbgp_attack CryptoAGG)
