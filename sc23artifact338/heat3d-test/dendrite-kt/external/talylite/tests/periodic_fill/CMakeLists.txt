cmake_minimum_required(VERSION 3.5.0 FATAL_ERROR)

project(test_periodic_fill)

add_executable(periodic_test
  # header files (not needed, but added so they show up in IDEs)
  include/pbc_test_input_data.h

  # source files to compile
  src/periodic_test.cpp
)

# add the ./include directory to the include search path so we can #include
# the header files we put in include/
target_include_directories(periodic_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

if (NOT TARGET talyfem)
  find_package(talyfem REQUIRED)
endif()

# this will also copy all the compiler definitions, libraries, etc. that TalyFEM
# was configured with
target_link_libraries(periodic_test talyfem)

# add tests
enable_testing()
add_test(NAME "test_periodic_fill"
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_test.py
)
set_property(TEST "test_periodic_fill"
  PROPERTY ENVIRONMENT
  EXECUTABLE_PATH=$<TARGET_FILE:periodic_test>
  TEST_DATA_DIR=${CMAKE_CURRENT_SOURCE_DIR}/tests
  TALYFEM_DIR=$<TARGET_PROPERTY:talyfem,SOURCE_DIR>
)
