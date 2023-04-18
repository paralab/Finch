cmake_minimum_required(VERSION 3.5.0 FATAL_ERROR)

project(test_mesquite)

add_executable(mesquite_test
  # header files (not needed, but added so they show up in IDEs)

  # source files to compile
  src/main.cpp
)

# add the ./include directory to the include search path so we can #include
# the header files we put in include/
target_include_directories(mesquite_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

if (NOT TARGET talyfem)
  find_package(talyfem REQUIRED)
endif()

# this will also copy all the compiler definitions, libraries, etc. that TalyFEM
# was configured with
target_link_libraries(mesquite_test talyfem)

# add tests (no python, just test the executable runs)
enable_testing()
add_test(NAME "test_mesquite"
  COMMAND mesquite_test
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test
)
