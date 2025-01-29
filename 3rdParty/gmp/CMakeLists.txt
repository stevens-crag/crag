#Specify the version being used aswell as the language
cmake_minimum_required(VERSION 2.8.12)
cmake_policy(VERSION 2.8.12)

include(ExternalProject)

if(UNIX)
  ExternalProject_Add(m4
    URL http://ftp.gnu.org/gnu/m4/m4-1.4.19.tar.bz2
    URL_HASH MD5=a8ff6fc3235eab3a4adb53611fdddacd
    CONFIGURE_COMMAND "<SOURCE_DIR>/configure" "--prefix=<INSTALL_DIR>"
  )
  
  ExternalProject_Get_Property(m4 install_dir)
  ExternalProject_Add(gmp_external
    DEPENDS m4
    URL http://ftp.gnu.org/gnu/gmp/gmp-6.3.0.tar.bz2
    URL_HASH MD5=c1cd6ef33085e9cb818b9b08371f9000
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

  ExternalProject_Add(MPIR__EXTERNAL
    URL "https://github.com/BrianGladman/mpir/archive/4f6d27d18ab767f44f57641c51bc0fedc5d165dc.zip"
    URL_HASH "MD5=a469b0fc50de12eec779aa537797698d"
    CONFIGURE_COMMAND ""
    BUILD_IN_SOURCE 1
    BUILD_COMMAND "${CMAKE_COMMAND}" ARGS -E chdir "<SOURCE_DIR>/msvc/vs19" "${CMAKE_MAKE_PROGRAM}" "dll_mpir_gc\\dll_mpir_gc.vcxproj"  "/p:SolutionDir=<SOURCE_DIR>\\msvc\\vs19\\" "/p:configuration=$<CONFIG>" "/p:platform=x64" "/p:TargetName=mpir$<$<CONFIG:Debug>:d>"
    INSTALL_DIR "${GRD_EXTERNAL_INSTALL_DIR}"
    INSTALL_COMMAND ${CMAKE_COMMAND} -E copy_if_different <SOURCE_DIR>/dll/x64/$<CONFIG>/gmp.h <INSTALL_DIR>/include/gmp.h
    COMMAND ${CMAKE_COMMAND} -E copy_if_different <SOURCE_DIR>/dll/x64/$<CONFIG>/gmpxx.h <INSTALL_DIR>/include/gmpxx.h
    COMMAND ${CMAKE_COMMAND} -E copy_if_different <SOURCE_DIR>/dll/x64/$<CONFIG>/mpir.h <INSTALL_DIR>/include/mpir.h
    COMMAND ${CMAKE_COMMAND} -E copy <SOURCE_DIR>/msvc/vs19/dll_mpir_gc/x64/$<CONFIG>/mpir$<$<CONFIG:Debug>:d>.dll <INSTALL_DIR>/lib/mpir$<$<CONFIG:Debug>:d>.dll
    COMMAND ${CMAKE_COMMAND} -E copy <SOURCE_DIR>/msvc/vs19/dll_mpir_gc/x64/$<CONFIG>/mpir$<$<CONFIG:Debug>:d>.lib <INSTALL_DIR>/lib/mpir$<$<CONFIG:Debug>:d>.lib
    # needed because SOURCE_DIR is cleared when another config is built
    COMMAND ${CMAKE_COMMAND} -E copy_directory <SOURCE_DIR>/dll/x64/$<CONFIG>/ <INSTALL_DIR>/lib/mpir-build/$<CONFIG>/
    COMMAND ${CMAKE_COMMAND} -E copy <SOURCE_DIR>/msvc/vs19/dll_mpir_gc/x64/$<CONFIG>/mpir$<$<CONFIG:Debug>:d>.dll <INSTALL_DIR>/lib/mpir-build/$<CONFIG>/mpir$<$<CONFIG:Debug>:d>.dll
    COMMAND ${CMAKE_COMMAND} -E copy <SOURCE_DIR>/msvc/vs19/dll_mpir_gc/x64/$<CONFIG>/mpir$<$<CONFIG:Debug>:d>.lib <INSTALL_DIR>/lib/mpir-build/$<CONFIG>/mpir$<$<CONFIG:Debug>:d>.lib
  )

  set_target_properties(MPIR__EXTERNAL PROPERTIES FOLDER "3rdParty/gmp")

  ExternalProject_Get_Property(MPIR__EXTERNAL install_dir)

  file(MAKE_DIRECTORY "${install_dir}/include")

  add_library(gmp STATIC IMPORTED GLOBAL)
  add_dependencies(gmp MPIR__EXTERNAL)  
  set_target_properties(gmp PROPERTIES
    IMPORTED_CONFIGURATIONS "Debug;Release"
    IMPORTED_LOCATION_DEBUG "${install_dir}/lib/mpird.lib"
    IMPORTED_LOCATION_RELEASE "${install_dir}/lib/mpir.lib"
    INTERFACE_INCLUDE_DIRECTORIES "${install_dir}/include"
    IMPORTED_LINK_INTERFACE_LANGUAGES "C"
  )

  add_library(gmpxx STATIC IMPORTED GLOBAL)
  add_dependencies(gmpxx MPIR__EXTERNAL )
  set_target_properties(gmpxx PROPERTIES
    IMPORTED_CONFIGURATIONS "Debug;Release"
    IMPORTED_LOCATION_DEBUG "${install_dir}/lib/mpird.lib"
    IMPORTED_LOCATION_RELEASE "${install_dir}/lib/mpir.lib"
    INTERFACE_INCLUDE_DIRECTORIES "${install_dir}/include"
    IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
  )
  
  set(MPIR_INSTALL_DIR ${install_dir}/lib PARENT_SCOPE)
ELSE()
   message( FATAL_ERROR "You should use UNIX or MSVC" )
ENDIF()
