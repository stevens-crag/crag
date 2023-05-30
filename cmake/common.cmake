function(crag_library name)
  set(srcs)
  set(current_arg 1)
  while(current_arg LESS ${ARGC})
    list(APPEND srcs "src/${ARGV${current_arg}}.cpp")
    math(EXPR current_arg "${current_arg}+1")
  endwhile()
  file(GLOB headers "${CMAKE_CURRENT_SOURCE_DIR}/include/*.h")
  add_library(${name} STATIC ${srcs} ${headers})
  target_include_directories(${name} PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
  if (UNIX)
    target_compile_options(${name} PRIVATE -Wall)
  endif()
  set_target_properties(${name} PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY lib
  )
endfunction()

function(crag_main name)
  get_filename_component(LIBNAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  file(GLOB headers "${CMAKE_CURRENT_SOURCE_DIR}/include/*.h")
  add_executable("${LIBNAME}_main_${name}" "main/${name}.cpp")

  set_target_properties("${LIBNAME}_main_${name}" PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY bin
    OUTPUT_NAME ${name}
  )

  list(REMOVE_AT ARGV 0)
  foreach(lib ${ARGV})
    target_link_libraries("${LIBNAME}_main_${name}" PUBLIC ${lib})
  endforeach()

  if (UNIX)
    target_compile_options("${LIBNAME}_main_${name}" PRIVATE -Wall)
  endif()

endfunction()

function(crag_exe_path_with_dll_env target return_var)
  if(MSVC)
    set("${return_var}" "${CMAKE_COMMAND}" -E env
        "PATH=%PATH%$<SEMICOLON>${MPIR_INSTALL_DIR}"
        "$<TARGET_FILE:${target}>"
        PARENT_SCOPE)
   else()
    set("${return_var}" "$<TARGET_FILE:${target}>" PARENT_SCOPE)
  endif()
endfunction()


function(crag_test name)
  get_filename_component(LIBNAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)

  add_executable("${LIBNAME}_test_${name}" "test/${name}.cpp")
  target_link_libraries("${LIBNAME}_test_${name}" PRIVATE GTest::gtest_main)

  set_target_properties("${LIBNAME}_test_${name}" PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY test
    OUTPUT_NAME ${name}
    VS_DEBUGGER_ENVIRONMENT "PATH=${MPIR_INSTALL_DIR}"
  )

  list(REMOVE_AT ARGV 0)
  foreach(lib ${ARGV})
    target_link_libraries("${LIBNAME}_test_${name}" PUBLIC ${lib})
  endforeach()

  if (UNIX)
    target_compile_options("${LIBNAME}_test_${name}" PRIVATE -Wall)
  endif()

  crag_exe_path_with_dll_env("${LIBNAME}_test_${name}" test_path)
  add_test(NAME ${name} COMMAND ${test_path})

endfunction()
