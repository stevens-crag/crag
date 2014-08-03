#Specify the version being used aswell as the language
cmake_minimum_required(VERSION 2.8.12)
cmake_policy(VERSION 2.8.12)

include("../cmake/common.cmake")

crag_library(SLPv1
  WhiteheadGraph
  FreeGroup
  StraightLineProgramWord
  StraightLineProgramWord_operators
  StraightLineProgramWord_truncate
  StraightLineProgramWord_split
)

target_link_libraries(SLPv1
  PUBLIC Elt
  PRIVATE gmpxx
  PUBLIC Maps
  PRIVATE SbgpFG
)

crag_main(test_cs SLPv1)

crag_library(SLPv2
  EndomorphismSLP
  slp
  slp_recompression
)

target_link_libraries(SLPv2
  PUBLIC boost_pool
  PUBLIC crag_general
  PUBLIC gmpxx
  PRIVATE gmp
)

find_library(tcmalloc_minimal libtcmalloc_minimal.so)
if (tcmalloc_minimal)
  target_link_libraries(SLPv2 PRIVATE "${tcmalloc_minimal}")
endif()


crag_main(cache_length_expirement SLPv2 boost_pool_gmp_allocator)
crag_main(hash_reduce             SLPv2 boost_pool_gmp_allocator)
crag_main(profile_matching_new    SLPv2 boost_pool_gmp_allocator)
crag_main(profile_reduce          SLPv2 boost_pool_gmp_allocator)
crag_main(profile_reduce_narrow   SLPv2 boost_pool_gmp_allocator)
crag_main(reduce_structure        SLPv2 boost_pool_gmp_allocator)

crag_test(AAGCryptoTest           SLPv2)
crag_test(arithmetic_sequence     SLPv2)
crag_test(EndomorphismSLP_test    SLPv2)
crag_test(FGACryptoTest           SLPv2)
crag_test(permutation16           SLPv2)
crag_test(slp_common_prefix       SLPv2)
crag_test(slp_inspector           SLPv2)
crag_test(slp_pattern_matching    SLPv2)
crag_test(slp_recompression_test  SLPv2)
crag_test(slp_reduce              SLPv2)
crag_test(slp_vertex              SLPv2)
crag_test(slp_vertex_hash         SLPv2)
crag_test(slp_vertex_word         SLPv2)








#target_link_libraries(CryptoTripleDecomposition
#  PUBLIC BraidGroup
#  PRIVATE ranlib
#)
