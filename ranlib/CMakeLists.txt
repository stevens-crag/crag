cmake_minimum_required(VERSION 2.8.12)
cmake_policy(VERSION 2.8.12)

crag_library(ranlib
  RanlibCPP
  com
  dcdflib
  ipmpar
  linpack
  ranlib
)

crag_main(test_RandLib ranlib)
