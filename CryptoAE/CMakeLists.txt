#Specify the version being used aswell as the language
cmake_minimum_required(VERSION 2.8.12)
cmake_policy(VERSION 2.8.12)

include("../cmake/common.cmake")

crag_library(CryptoAE
  AEProtocol
  TTPAttack
)

target_link_libraries(CryptoAE
  PUBLIC BraidGroup
  PRIVATE ranlib
)

crag_main(test_TTP CryptoAE ranlib Boost::program_options)

crag_test(test_cryptoae CryptoAE)
