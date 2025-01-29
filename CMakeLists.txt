cmake_minimum_required(VERSION 3.8)

### Hunter snapshot that will be used ###
include("cmake/HunterGate.cmake")
HunterGate(
  URL "https://github.com/cpp-pm/hunter/archive/refs/heads/master.zip"
  SHA1 "90cb83da9a08f7ba0acf5028f6c48deaffff1760"
)

project(crag)

### Download dependencies ###
hunter_add_package(Boost COMPONENTS program_options)
hunter_add_package(GTest)
hunter_add_package(benchmark)

### Find dependencies ###
find_package(GTest CONFIG REQUIRED) # GTest
find_package(Boost CONFIG REQUIRED program_options) # Boost::boost, Boost::program_options
find_package(benchmark CONFIG REQUIRED) # benchmark::benchmark

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

set_property(GLOBAL PROPERTY USE_FOLDERS ON) #turn on folders for VS

# Add subdirectories (All new directories should go here)
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
add_subdirectory(FiniteField)
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
add_subdirectory(Polynomial)
add_subdirectory(Matrix)
add_subdirectory(Random)
add_subdirectory(Kayawood)
add_subdirectory(Walnut)
