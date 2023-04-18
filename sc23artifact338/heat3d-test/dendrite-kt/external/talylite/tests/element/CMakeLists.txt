cmake_minimum_required(VERSION 3.5.0 FATAL_ERROR)

project(test_element)

add_executable(element_test
  # header files (not needed, but added so they show up in IDEs)
  include/element_test.h
  include/globals.h

  include/tests/elem_1d.h
  include/tests/elem_2d_box.h
  include/tests/elem_2d_triangle.h
  include/tests/elem_3d_hex.h
  include/tests/elem_3d_tet.h

  # source files to compile
  src/main.cpp
  src/element_test.cpp
)

# add the ./include directory to the include search path so we can #include
# the header files we put in include/
target_include_directories(element_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

if (NOT TARGET talyfem)
  find_package(talyfem REQUIRED)
endif()

# this will also copy all the compiler definitions, libraries, etc. that TalyFEM
# was configured with
target_link_libraries(element_test talyfem)

# add tests (no python, just test the executable runs)
enable_testing()
add_test(NAME "test_element"
  COMMAND element_test
)
