#Specify the version being used aswell as the language
cmake_minimum_required(VERSION 2.8.12)
cmake_policy(VERSION 2.8.12)

include(ExternalProject)

if(UNIX)
  ExternalProject_Add(m4
    URL http://ftp.gnu.org/gnu/m4/m4-1.4.17.tar.bz2 
    URL_HASH MD5=8a1787edcba75ae5cd1dc40d7d8ed03a
    CONFIGURE_COMMAND "<SOURCE_DIR>/configure" "--prefix=<INSTALL_DIR>"
  )
  
  ExternalProject_Get_Property(m4 install_dir)
  ExternalProject_Add(gmp_external
    DEPENDS m4
    URL http://ftp.gnu.org/gnu/gmp/gmp-6.0.0a.tar.bz2
    URL_HASH MD5=b7ff2d88cae7f8085bd5006096eed470
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
  ExternalProject_Add(mpir-debug
    PREFIX mpir
    URL http://www.mpir.org/mpir-2.6.0.tar.bz2
    URL_HASH MD5=81be2f605330d030956c1ff515c58564
    SOURCE_DIR "src/mpir"
    PATCH_COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/mpir_zero_patched.c <SOURCE_DIR>/mpn/generic/zero.c
	  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/mpirxx_patched.h <SOURCE_DIR>/mpirxx.h
    CONFIGURE_COMMAND ""
    BUILD_IN_SOURCE 1
    BUILD_COMMAND ${CMAKE_COMMAND} -E chdir "<SOURCE_DIR>/build.vc10" ${CMAKE_MAKE_PROGRAM} lib_mpir_gc\\lib_mpir_gc.vcxproj /p:configuration=Debug /p:platform=x64 /p:SolutionDir=<SOURCE_DIR>\\build.vc10\\\\ /p:PlatformToolset=v120
    INSTALL_DIR "mpir"
    INSTALL_COMMAND ${CMAKE_COMMAND} -E copy "<SOURCE_DIR>/lib/x64/Debug/gmp.h" "<INSTALL_DIR>/include/gmp.h"
      COMMAND ${CMAKE_COMMAND} -E copy_directory "<SOURCE_DIR>/lib/x64/Debug" "<INSTALL_DIR>/lib/x64/Debug"
      COMMAND ${CMAKE_COMMAND} -E copy "<SOURCE_DIR>/lib/x64/Debug/mpir.lib" "<INSTALL_DIR>/lib/mpir-debug.lib"
      COMMAND ${CMAKE_COMMAND} -E copy "<SOURCE_DIR>/lib/x64/Debug/mpir.pdb" "<INSTALL_DIR>/lib/mpir-debug.pdb"
  )
  
  add_library(gmp STATIC IMPORTED GLOBAL)
  add_library(gmpxx STATIC IMPORTED GLOBAL)
  ExternalProject_Get_Property(mpir-debug install_dir)

  file(MAKE_DIRECTORY "${install_dir}/include")
  
  set_target_properties(gmp PROPERTIES
    IMPORTED_CONFIGURATIONS "Debug;Release"
    IMPORTED_LOCATION_DEBUG "${install_dir}/lib/mpir-debug.lib"
    INTERFACE_INCLUDE_DIRECTORIES "${install_dir}/include"
  )
  
  add_dependencies(gmp mpir-debug)
  
  ExternalProject_Get_Property(mpir-debug source_dir)
  
  ExternalProject_Add(mpir-release
    PREFIX mpir
    DEPENDS mpir-debug
    DOWNLOAD_COMMAND ""
    SOURCE_DIR "src/mpir"
    CONFIGURE_COMMAND ""
    BUILD_IN_SOURCE 1
    BUILD_COMMAND ${CMAKE_COMMAND} -E chdir "<SOURCE_DIR>/build.vc10" ${CMAKE_MAKE_PROGRAM} lib_mpir_gc\\lib_mpir_gc.vcxproj /p:configuration=Release /p:platform=x64 /p:SolutionDir=<SOURCE_DIR>\\build.vc10\\\\ /p:PlatformToolset=v120
    INSTALL_DIR "mpir"
    INSTALL_COMMAND ${CMAKE_COMMAND} -E copy "<SOURCE_DIR>/lib/x64/Release/gmp.h" "<INSTALL_DIR>/include/gmp.h"
	  COMMAND ${CMAKE_COMMAND} -E copy "<SOURCE_DIR>/lib/x64/Release/mpir.h" "<INSTALL_DIR>/include/mpir.h"
      COMMAND ${CMAKE_COMMAND} -E copy "<SOURCE_DIR>/lib/x64/Release/gmpxx.h" "<INSTALL_DIR>/include/gmpxx.h"
      COMMAND ${CMAKE_COMMAND} -E copy_directory "<SOURCE_DIR>/lib/x64/Release" "<INSTALL_DIR>/lib/x64/Release"
      COMMAND ${CMAKE_COMMAND} -E copy "<SOURCE_DIR>/lib/x64/Release/mpir.lib" "<INSTALL_DIR>/lib/mpir-release.lib"
      COMMAND ${CMAKE_COMMAND} -E copy "<SOURCE_DIR>/lib/x64/Release/mpir.pdb" "<INSTALL_DIR>/lib/mpir-release.pdb"
  )
  
  ExternalProject_Get_Property(mpir-release install_dir)
  set_property(TARGET gmp 
    PROPERTY IMPORTED_LOCATION_RELEASE "${install_dir}/lib/mpir-release.lib"
  )
  
  add_dependencies(gmp mpir-release)

  ExternalProject_Add(mpirxx-debug
    PREFIX mpir
    DEPENDS mpir-debug
    DOWNLOAD_COMMAND ""
    SOURCE_DIR "src/mpir"
    CONFIGURE_COMMAND ""
    BUILD_IN_SOURCE 1
    BUILD_COMMAND ${CMAKE_COMMAND} -E chdir "<SOURCE_DIR>/build.vc10" ${CMAKE_MAKE_PROGRAM} lib_mpir_cxx\\lib_mpir_cxx.vcxproj /p:configuration=Debug /p:platform=x64 /p:SolutionDir=<SOURCE_DIR>\\build.vc10\\\\ /p:PlatformToolset=v120
    INSTALL_DIR "mpir"
    INSTALL_COMMAND ${CMAKE_COMMAND} -E copy "<SOURCE_DIR>/lib/x64/Debug/gmpxx.h" "<INSTALL_DIR>/include/gmpxx.h"
      COMMAND ${CMAKE_COMMAND} -E copy_directory "<SOURCE_DIR>/lib/x64/Debug" "<INSTALL_DIR>/lib/x64/Debug"
      COMMAND ${CMAKE_COMMAND} -E copy "<SOURCE_DIR>/lib/x64/Debug/mpirxx.lib" "<INSTALL_DIR>/lib/mpirxx-debug.lib"
      COMMAND ${CMAKE_COMMAND} -E copy "<SOURCE_DIR>/lib/x64/Debug/mpirxx.pdb" "<INSTALL_DIR>/lib/mpirxx-debug.pdb"
  )

  ExternalProject_Add(mpirxx-release
    PREFIX mpir
    DEPENDS mpir-debug
    DOWNLOAD_COMMAND ""
    SOURCE_DIR "src/mpir"
    CONFIGURE_COMMAND ""
    BUILD_IN_SOURCE 1
    BUILD_COMMAND ${CMAKE_COMMAND} -E chdir "<SOURCE_DIR>/build.vc10" ${CMAKE_MAKE_PROGRAM} lib_mpir_cxx\\lib_mpir_cxx.vcxproj /p:configuration=Release /p:platform=x64 /p:SolutionDir=<SOURCE_DIR>\\build.vc10\\\\ /p:PlatformToolset=v120
    INSTALL_DIR "mpir"
    INSTALL_COMMAND ${CMAKE_COMMAND} -E copy "<SOURCE_DIR>/lib/x64/Release/gmpxx.h" "<INSTALL_DIR>/include/gmpxx.h"
      COMMAND ${CMAKE_COMMAND} -E copy_directory "<SOURCE_DIR>/lib/x64/Release" "<INSTALL_DIR>/lib/x64/Release"
      COMMAND ${CMAKE_COMMAND} -E copy "<SOURCE_DIR>/lib/x64/Release/mpirxx.lib" "<INSTALL_DIR>/lib/mpirxx-release.lib"
      COMMAND ${CMAKE_COMMAND} -E copy "<SOURCE_DIR>/lib/x64/Release/mpirxx.pdb" "<INSTALL_DIR>/lib/mpirxx-release.pdb"
  )


  set_target_properties(gmpxx PROPERTIES
    IMPORTED_CONFIGURATIONS "Debug;Release"
    INTERFACE_INCLUDE_DIRECTORIES "${install_dir}/include"
    IMPORTED_LOCATION_DEBUG "${install_dir}/lib/mpirxx-debug.lib"
    IMPORTED_LOCATION_RELEASE "${install_dir}/lib/mpirxx-release.lib"
  )

  add_dependencies(gmpxx mpirxx-debug mpirxx-release)

    
ELSE()
   message( FATAL_ERROR "You should use UNIX or MSVC" )
ENDIF()
