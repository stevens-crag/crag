cmake_minimum_required(VERSION 3.8)

### Hunter snapshot that will be used ###
include("cmake/HunterGate.cmake")
HunterGate(
  URL "https://github.com/ruslo/hunter/archive/v0.20.9.tar.gz"
  SHA1 "8f7cf6d719d05c49224a4ebeba90c6bb90b90f4a"
)

project(crag)

### Download dependencies ###
hunter_add_package(Boost)
hunter_add_package(GTest)

### Find dependencies ###
find_package(GTest CONFIG REQUIRED) # GTest::main
find_package(Boost CONFIG REQUIRED) # Boost::boost, header only libs

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")

IF(UNIX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")

  set(CMAKE_CXX_FLAGS_ASAN "-g -O3 -fsanitize=address,undefined"
          CACHE STRING "Compiler flags in asan build"
          FORCE)

  set(CMAKE_CXX_FLAGS_TSAN "-g -O3 -fsanitize=thread"
          CACHE STRING "Compiler flags in tsan build"
          FORCE)

ENDIF(UNIX)

enable_testing()

if(BUILD_SHARED_LIBS)
  message(FATAL_ERROR "Shared build is not supported")
endif()

IF(MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W3 -EHsc -nologo")
  foreach (flag_var
    CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE
    CMAKE_CXX_FLAGS_MINSIZEREL CMAKE_CXX_FLAGS_RELWITHDEBINFO)
    if(NOT BUILD_SHARED_LIBS)
      #force using static runtime
      string(REPLACE "/MD" "-MT" ${flag_var} "${${flag_var}}")
    endif()
  endforeach()
ENDIF(MSVC)

set_property(GLOBAL PROPERTY USE_FOLDERS ON) #turn on folders for VS

# Add subdirectories (All new directories should go here)
add_subdirectory(3rdParty/argagg)
add_subdirectory(3rdParty/gmp)
add_subdirectory(Alphabet)
add_subdirectory(BraidGroup)
add_subdirectory(CryptoAAG)
add_subdirectory(CryptoAE)
add_subdirectory(CryptoKL)
add_subdirectory(CryptoShftConj)
add_subdirectory(CryptoTripleDecomposition)
add_subdirectory(Elt)
add_subdirectory(Equation)
add_subdirectory(Examples)
add_subdirectory(FreeGroup)
add_subdirectory(general)
add_subdirectory(Graph)
add_subdirectory(Graphics)
add_subdirectory(Group)
add_subdirectory(HigmanGroup)
add_subdirectory(Maps)
add_subdirectory(ranlib)
add_subdirectory(StringSimilarity)
add_subdirectory(SbgpFG)
add_subdirectory(TheGrigorchukGroup)
