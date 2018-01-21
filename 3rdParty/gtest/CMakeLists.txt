cmake_minimum_required(VERSION 2.8.12)
cmake_policy(VERSION 2.8.12)

find_package(Threads REQUIRED)
include(ExternalProject)

if(MSVC)
  set(SILENCE_DEPRECATION "-DCMAKE_CXX_FLAGS=/D_SILENCE_TR1_NAMESPACE_DEPRECATION_WARNING")
endif()

ExternalProject_Add(
    gtest_build
    URL https://github.com/google/googletest/archive/release-1.8.0.zip 
    URL_HASH MD5=adfafc8512ab65fd3cf7955ef0100ff5
    CMAKE_ARGS -DBUILD_GTEST=1 -DBUILD_GMOCK=0 -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR> -DCMAKE_DEBUG_POSTFIX:STRING=d ${SILENCE_DEPRECATION}
)

ExternalProject_Get_Property(gtest_build install_dir)

add_library(gtest STATIC IMPORTED GLOBAL)

#workaround for http://public.kitware.com/Bug/bug_relationship_graph.php?bug_id=14495
file(MAKE_DIRECTORY "${install_dir}/include")

if(UNIX)
	set_target_properties(gtest PROPERTIES
	  IMPORTED_LOCATION "${install_dir}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}gtest${CMAKE_STATIC_LIBRARY_SUFFIX}"
	  INTERFACE_INCLUDE_DIRECTORIES "${install_dir}/include"
	  INTERFACE_LINK_LIBRARIES "${CMAKE_THREAD_LIBS_INIT}"
	)
else()
	set_target_properties(gtest PROPERTIES
	  IMPORTED_LOCATION_DEBUG "${install_dir}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}gtestd${CMAKE_STATIC_LIBRARY_SUFFIX}"
	  IMPORTED_LOCATION_RELEASE "${install_dir}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}gtest${CMAKE_STATIC_LIBRARY_SUFFIX}"
	  INTERFACE_INCLUDE_DIRECTORIES "${install_dir}/include"
	  INTERFACE_LINK_LIBRARIES "${CMAKE_THREAD_LIBS_INIT}"
	)
endif()

add_library(gtest_main STATIC IMPORTED GLOBAL)
set_target_properties(gtest_main PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "${install_dir}/include"
  INTERFACE_LINK_LIBRARIES "${CMAKE_THREAD_LIBS_INIT};gtest"
)

if(UNIX)
	set_target_properties(gtest_main PROPERTIES
      IMPORTED_LOCATION "${install_dir}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}gtest_main${CMAKE_STATIC_LIBRARY_SUFFIX}"
	)
else()
	set_target_properties(gtest_main PROPERTIES
	  IMPORTED_LOCATION_DEBUG "${install_dir}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}gtest_maind${CMAKE_STATIC_LIBRARY_SUFFIX}"
	  IMPORTED_LOCATION_RELEASE "${install_dir}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}gtest_main${CMAKE_STATIC_LIBRARY_SUFFIX}"
	)
endif()

add_dependencies(gtest gtest_build)
add_dependencies(gtest_main gtest_build)
