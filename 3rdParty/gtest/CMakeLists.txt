cmake_minimum_required(VERSION 3.0.0)
cmake_policy(VERSION 3.0.0)

find_package(Threads REQUIRED)
include(ExternalProject)

ExternalProject_Add(
    gtest_build
    URL http://googletest.googlecode.com/files/gtest-1.7.0.zip
    URL_HASH MD5=2d6ec8ccdf5c46b05ba54a9fd1d130d7
    INSTALL_COMMAND ""
)

ExternalProject_Get_Property(gtest_build source_dir)
ExternalProject_Get_Property(gtest_build binary_dir)

add_library(gtest STATIC IMPORTED GLOBAL)

#workaround for http://public.kitware.com/Bug/bug_relationship_graph.php?bug_id=14495
file(MAKE_DIRECTORY "${source_dir}/include")

set_target_properties(gtest PROPERTIES
  IMPORTED_LOCATION "${binary_dir}/${CMAKE_STATIC_LIBRARY_PREFIX}gtest${CMAKE_STATIC_LIBRARY_SUFFIX}"
  INTERFACE_INCLUDE_DIRECTORIES "${source_dir}/include"
  INTERFACE_LINK_LIBRARIES "${CMAKE_THREAD_LIBS_INIT}"
)


add_library(gtest_main STATIC IMPORTED GLOBAL)
set_target_properties(gtest_main PROPERTIES
  IMPORTED_LOCATION "${binary_dir}/${CMAKE_STATIC_LIBRARY_PREFIX}gtest_main${CMAKE_STATIC_LIBRARY_SUFFIX}"
  INTERFACE_INCLUDE_DIRECTORIES "${source_dir}/include"
  INTERFACE_LINK_LIBRARIES "${CMAKE_THREAD_LIBS_INIT};gtest"
)

add_dependencies(gtest gtest_build)
add_dependencies(gtest_main gtest_build)
