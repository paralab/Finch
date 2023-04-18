cmake_minimum_required(VERSION 3.5.0 FATAL_ERROR)

project(test_pbc_dd_map)

add_executable(pbc_dd_index_test
  # header files (not needed, but added so they show up in IDEs)
  include/pbc_equation.h
  include/pbc_input_data.h

  # source files to compile
  src/pbc_dd_index_test.cpp
)

# add the ./include directory to the include search path so we can #include
# the header files we put in include/
target_include_directories(pbc_dd_index_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

if (NOT TARGET talyfem)
  find_package(talyfem REQUIRED)
endif()

# this will also copy all the compiler definitions, libraries, etc. that TalyFEM
# was configured with
target_link_libraries(pbc_dd_index_test talyfem)

# add tests
enable_testing()
add_test(NAME "test_pbc_dd_map"
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_test.py
)
set_property(TEST "test_pbc_dd_map"
  PROPERTY ENVIRONMENT
  EXECUTABLE_PATH=$<TARGET_FILE:pbc_dd_index_test>
  TEST_DATA_DIR=${CMAKE_CURRENT_SOURCE_DIR}/tests
  TALYFEM_DIR=$<TARGET_PROPERTY:talyfem,SOURCE_DIR>
)
