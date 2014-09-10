cmake_minimum_required(VERSION 2.8.12)
cmake_policy(VERSION 2.8.12)

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

if(UNIX)
	set_target_properties(gtest PROPERTIES
	  IMPORTED_LOCATION "${binary_dir}/${CMAKE_STATIC_LIBRARY_PREFIX}gtest${CMAKE_STATIC_LIBRARY_SUFFIX}"
	  INTERFACE_INCLUDE_DIRECTORIES "${source_dir}/include"
	  INTERFACE_LINK_LIBRARIES "${CMAKE_THREAD_LIBS_INIT}"
	)
else()
	set_target_properties(gtest PROPERTIES
	  IMPORTED_LOCATION_DEBUG "${binary_dir}/Debug/${CMAKE_STATIC_LIBRARY_PREFIX}gtest${CMAKE_STATIC_LIBRARY_SUFFIX}"
	  IMPORTED_LOCATION_RELEASE "${binary_dir}/Release/${CMAKE_STATIC_LIBRARY_PREFIX}gtest${CMAKE_STATIC_LIBRARY_SUFFIX}"
	  INTERFACE_INCLUDE_DIRECTORIES "${source_dir}/include"
	  INTERFACE_LINK_LIBRARIES "${CMAKE_THREAD_LIBS_INIT}"
	)
endif()



add_library(gtest_main STATIC IMPORTED GLOBAL)
set_target_properties(gtest_main PROPERTIES
  IMPORTED_LOCATION "${binary_dir}/${CMAKE_STATIC_LIBRARY_PREFIX}gtest_main${CMAKE_STATIC_LIBRARY_SUFFIX}"
  INTERFACE_INCLUDE_DIRECTORIES "${source_dir}/include"
  INTERFACE_LINK_LIBRARIES "${CMAKE_THREAD_LIBS_INIT};gtest"
)

if(UNIX)
	set_target_properties(gtest_main PROPERTIES
      IMPORTED_LOCATION "${binary_dir}/${CMAKE_STATIC_LIBRARY_PREFIX}gtest_main${CMAKE_STATIC_LIBRARY_SUFFIX}"
      INTERFACE_INCLUDE_DIRECTORIES "${source_dir}/include"
      INTERFACE_LINK_LIBRARIES "${CMAKE_THREAD_LIBS_INIT};gtest"
	)
else()
	set_target_properties(gtest_main PROPERTIES
	  IMPORTED_LOCATION_DEBUG "${binary_dir}/Debug/${CMAKE_STATIC_LIBRARY_PREFIX}gtest_main${CMAKE_STATIC_LIBRARY_SUFFIX}"
	  IMPORTED_LOCATION_RELEASE "${binary_dir}/Release/${CMAKE_STATIC_LIBRARY_PREFIX}gtest_main${CMAKE_STATIC_LIBRARY_SUFFIX}"
	  INTERFACE_INCLUDE_DIRECTORIES "${source_dir}/include"
	  INTERFACE_LINK_LIBRARIES "${CMAKE_THREAD_LIBS_INIT}"
	)
endif()

add_dependencies(gtest gtest_build)
add_dependencies(gtest_main gtest_build)
