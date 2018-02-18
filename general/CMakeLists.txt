cmake_minimum_required(VERSION 3.8)

include(../cmake/common.cmake)

crag_library(boost_pool_gmp_allocator
  gmp_boost_pool_allocator
)

target_link_libraries(boost_pool_gmp_allocator
  PRIVATE Boost::boost
  PRIVATE gmp
)

crag_library(crag_general
  Permutation
  ConfigFile
  BalancedTree
  VectorEnumerator
)

target_link_libraries(crag_general PRIVATE ranlib)

crag_main(mask crag_general ranlib)
