cmake_minimum_required(VERSION 2.8.12)
cmake_policy(VERSION 2.8.12)

include(ExternalProject)

ExternalProject_Add(
    argagg_build
    URL https://github.com/vietjtnguyen/argagg/archive/0.4.6.zip
    URL_HASH MD5=221e5df3d47a1305a347cea3f0ae9fed
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR> -DCMAKE_DEBUG_POSTFIX:STRING=d
)

ExternalProject_Get_Property(argagg_build install_dir)

add_library(argagg INTERFACE IMPORTED GLOBAL)

#workaround for http://public.kitware.com/Bug/bug_relationship_graph.php?bug_id=14495
file(MAKE_DIRECTORY "${install_dir}/include")

set_target_properties(argagg PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "${install_dir}/include"
)

add_dependencies(argagg argagg_build)
