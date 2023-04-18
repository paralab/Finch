cmake_minimum_required(VERSION 3.5.0 FATAL_ERROR)

project(test_integrator)

add_executable(integrator_test
  # header files (not needed, but added so they show up in IDEs)
  include/integrator_test_grid_field.h
  include/integrator_test_input_data.h
  include/integrator_test_node_data.h
  include/volume_function.h

  # source files to compile
  src/main.cpp
  src/volume_function.cpp
)

# add the ./include directory to the include search path so we can #include
# the header files we put in include/
target_include_directories(integrator_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

if (NOT TARGET talyfem)
  find_package(talyfem REQUIRED)
endif()

# this will also copy all the compiler definitions, libraries, etc. that TalyFEM
# was configured with
target_link_libraries(integrator_test talyfem)

# add tests
enable_testing()
add_test(NAME "test_integrator"
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_test.py
)
set_property(TEST "test_integrator"
  PROPERTY ENVIRONMENT
  EXECUTABLE_PATH=$<TARGET_FILE:integrator_test>
  TEST_DATA_DIR=${CMAKE_CURRENT_SOURCE_DIR}/tests
  TALYFEM_DIR=$<TARGET_PROPERTY:talyfem,SOURCE_DIR>
)
