#Specify the version being used aswell as the language
cmake_minimum_required(VERSION 2.8.12)
cmake_policy(VERSION 2.8.12)

include(ExternalProject)

if(UNIX)
  ExternalProject_Add(m4
    URL http://ftp.gnu.org/gnu/m4/m4-1.4.18.tar.bz2 
    URL_HASH MD5=199fff9c87f1c920dab5c8757811e146
    PATCH_COMMAND patch -p0 -d "<SOURCE_DIR>" < "${CMAKE_CURRENT_LIST_DIR}/secure_snprintf.patch"
    CONFIGURE_COMMAND "<SOURCE_DIR>/configure" "--prefix=<INSTALL_DIR>"
  )
  
  ExternalProject_Get_Property(m4 install_dir)
  ExternalProject_Add(gmp_external
    DEPENDS m4
    URL http://ftp.gnu.org/gnu/gmp/gmp-6.1.2.tar.bz2
    URL_HASH MD5=8ddbb26dc3bd4e2302984debba1406a5
    CONFIGURE_COMMAND M4="${install_dir}/bin/m4" "<SOURCE_DIR>/configure" --enable-shared=off "--prefix=<INSTALL_DIR>" --enable-cxx
  )

  ExternalProject_Get_Property(gmp_external install_dir)
  file(MAKE_DIRECTORY "${install_dir}/include")
  
  add_library(gmp STATIC IMPORTED GLOBAL)
  add_dependencies(gmp gmp_external)
  
  set_target_properties(gmp PROPERTIES
    IMPORTED_LOCATION "${install_dir}/lib/libgmp.a"
    INTERFACE_INCLUDE_DIRECTORIES "${install_dir}/include"
  )
  
  add_library(gmpxx STATIC IMPORTED GLOBAL)
  add_dependencies(gmpxx gmp_external)
  
  set_target_properties(gmpxx PROPERTIES
    IMPORTED_LOCATION "${install_dir}/lib/libgmpxx.a"
    INTERFACE_INCLUDE_DIRECTORIES "${install_dir}/include"
  )
elseif(MSVC)
  if(NOT CMAKE_SIZEOF_VOID_P EQUAL 8 )
    message( FATAL_ERROR "Please, build with 64-bit VS generator" )
  endif()

  ExternalProject_Add(_mpir_external
    PREFIX mpir
    URL http://mpir.org/mpir-2.7.0.tar.bz2
    URL_HASH MD5=4aec04108dd24c49d45beb959d837b40
    SOURCE_DIR "src/mpir"
    CONFIGURE_COMMAND ""
    BUILD_IN_SOURCE 1
    BUILD_COMMAND "${CMAKE_COMMAND}" ARGS -E chdir "<SOURCE_DIR>/build.vc14" "${CMAKE_MAKE_PROGRAM}" "lib_mpir_gc\\lib_mpir_gc.vcxproj"  "/p:SolutionDir=<SOURCE_DIR>\\build.vc14\\x" "/p:configuration=Debug" "/p:platform=x64"
	  COMMAND "${CMAKE_COMMAND}" ARGS -E chdir "<SOURCE_DIR>/build.vc14" "${CMAKE_MAKE_PROGRAM}" "lib_mpir_gc\\lib_mpir_gc.vcxproj"  "/p:SolutionDir=<SOURCE_DIR>\\build.vc14\\x" "/p:configuration=Release" "/p:platform=x64"
	  COMMAND "${CMAKE_COMMAND}" ARGS -E chdir "<SOURCE_DIR>/build.vc14" "${CMAKE_MAKE_PROGRAM}" "lib_mpir_cxx\\lib_mpir_cxx.vcxproj"  "/p:SolutionDir=<SOURCE_DIR>\\build.vc14\\x" "/p:configuration=Debug" "/p:platform=x64"
	  COMMAND "${CMAKE_COMMAND}" ARGS -E chdir "<SOURCE_DIR>/build.vc14" "${CMAKE_MAKE_PROGRAM}" "lib_mpir_cxx\\lib_mpir_cxx.vcxproj"  "/p:SolutionDir=<SOURCE_DIR>\\build.vc14\\x" "/p:configuration=Release" "/p:platform=x64"
	INSTALL_DIR "mpir"
    INSTALL_COMMAND ${CMAKE_COMMAND} -E copy <SOURCE_DIR>/lib/xx64/Debug/mpir.h <INSTALL_DIR>/include/mpir.h
	  COMMAND ${CMAKE_COMMAND} -E copy <SOURCE_DIR>/lib/xx64/Debug/gmp.h <INSTALL_DIR>/include/gmp.h
      COMMAND ${CMAKE_COMMAND} -E copy <SOURCE_DIR>/lib/xx64/Debug/gmpxx.h <INSTALL_DIR>/include/gmpxx.h
      COMMAND ${CMAKE_COMMAND} -E copy_directory <SOURCE_DIR>/lib/xx64/Debug <INSTALL_DIR>/lib/x64/Debug
      COMMAND ${CMAKE_COMMAND} -E copy_directory <SOURCE_DIR>/lib/xx64/Release <INSTALL_DIR>/lib/x64/Release
  )

  add_library(gmp STATIC IMPORTED GLOBAL)
  add_dependencies(gmp _mpir_external)

  ExternalProject_Get_Property(_mpir_external install_dir)

  file(MAKE_DIRECTORY "${install_dir}/include")
  
  set_target_properties(gmp PROPERTIES
    IMPORTED_CONFIGURATIONS "Debug;Release"
    IMPORTED_LOCATION_DEBUG "${install_dir}/lib/x64/Debug/mpir.lib"
	IMPORTED_LOCATION_RELEASE "${install_dir}/lib/x64/Release/mpir.lib"
    INTERFACE_INCLUDE_DIRECTORIES "${install_dir}/include"
	IMPORTED_LINK_INTERFACE_LANGUAGES "C"
  )

  add_library(gmpxx STATIC IMPORTED GLOBAL)
  add_dependencies(gmpxx _mpir_external)
  set_target_properties(gmpxx PROPERTIES
    IMPORTED_CONFIGURATIONS "Debug;Release"
    IMPORTED_LOCATION_DEBUG "${install_dir}/lib/x64/Debug/mpirxx.lib"
	IMPORTED_LOCATION_RELEASE "${install_dir}/lib/x64/Release/mpirxx.lib"
    INTERFACE_INCLUDE_DIRECTORIES "${install_dir}/include"
	IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
  )
ELSE()
   message( FATAL_ERROR "You should use UNIX or MSVC" )
ENDIF()
