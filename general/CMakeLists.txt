#Specify the version being used aswell as the language
cmake_minimum_required(VERSION 2.8.12)
cmake_policy(VERSION 2.8.12)

include(ExternalProject)

ExternalProject_Add(boost_pool_external
  URL "${CMAKE_CURRENT_SOURCE_DIR}/boost_pool.tar.bz2"
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND ""
)
  
ExternalProject_Get_Property(boost_pool_external source_dir)
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/dumb.cpp" "namespace {void dumb(){}}")
add_library(boost_pool STATIC "${CMAKE_CURRENT_BINARY_DIR}/dumb.cpp")

target_include_directories(boost_pool INTERFACE ${source_dir})
target_compile_definitions(boost_pool INTERFACE BOOST_DISABLE_THREADS)
add_dependencies(boost_pool boost_pool_external)


include(../cmake/common.cmake)

crag_library(boost_pool_gmp_allocator
  gmp_boost_pool_allocator
)

target_link_libraries(boost_pool_gmp_allocator
  PRIVATE boost_pool
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
