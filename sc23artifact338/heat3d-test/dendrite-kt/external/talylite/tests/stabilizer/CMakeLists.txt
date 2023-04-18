cmake_minimum_required(VERSION 3.5.0 FATAL_ERROR)

project(test_stabilizer)

add_executable(stabilizer_test
  # header files (not needed, but added so they show up in IDEs)
  include/stabilizertest.h
  include/globals.h

  include/tests/box1d_cubic.h
  include/tests/box1d_linear.h
  include/tests/box1d_quadratic.h
  include/tests/box2d_cubic.h
  include/tests/box2d_linear.h
  include/tests/box2d_quadratic.h
  include/tests/box3d_cubic.h
  include/tests/box3d_linear.h
  include/tests/box3d_quadratic.h
  include/tests/tet3d_linear.h
  include/tests/tri2d_linear.h

  # source files to compile
  src/main.cpp
  src/stabilizertest.cpp
)

# add the ./include directory to the include search path so we can #include
# the header files we put in include/
target_include_directories(stabilizer_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

if (NOT TARGET talyfem)
  find_package(talyfem REQUIRED)
endif()

# this will also copy all the compiler definitions, libraries, etc. that TalyFEM
# was configured with
target_link_libraries(stabilizer_test talyfem)

# add tests (no python, just test the executable runs)
enable_testing()
add_test(NAME "test_stabilizer"
  COMMAND stabilizer_test
)
