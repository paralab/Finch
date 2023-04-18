cmake_minimum_required(VERSION 3.5.0 FATAL_ERROR)

project(test_basic_io)

add_executable(basic_io
  # header files (not needed, but added so they show up in IDEs)

  # source files to compile
  src/basic_io_test.cpp
)

# add the ./include directory to the include search path so we can #include
# the header files we put in include/
target_include_directories(basic_io PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

if (NOT TARGET talyfem)
  find_package(talyfem REQUIRED)
endif()

# this will also copy all the compiler definitions, libraries, etc. that TalyFEM
# was configured with
target_link_libraries(basic_io talyfem)

# add tests
enable_testing()
add_test(NAME "test_basic_io"
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_test.py
)
set_property(TEST "test_basic_io"
  PROPERTY ENVIRONMENT
  EXECUTABLE_PATH=$<TARGET_FILE:basic_io>
  TEST_DATA_DIR=${CMAKE_CURRENT_SOURCE_DIR}/tests
  TALYFEM_DIR=$<TARGET_PROPERTY:talyfem,SOURCE_DIR>
)
