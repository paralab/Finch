cmake_minimum_required(VERSION 3.5.0 FATAL_ERROR)

project(transient_heat)

add_executable(ht
  # header files (not needed, but added so they show up in IDEs)
  include/HTAnalyticSolution.h
  include/HTEquation.h
  include/HTGridField.h
  include/HTInputData.h
  include/HTNodeData.h

  # source files to compile
  src/main.cpp
)

# add the ./include directory to the include search path so we can #include
# the header files we put in include/
target_include_directories(ht PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

if (NOT TARGET talyfem)
  find_package(talyfem REQUIRED)
  include(TalyDoxygen)  # for add_doxygen()
  include(TalyReproducibility)  # for add_taly_reproducibility()
endif()

# this will also copy all the compiler definitions, etc. that TalyFEM was
# compiled with (that were marked as PUBLIC) - e.g. this will automatically
# add "-DENABLE_HDF5" if TalyFEM was compiled with HDF5 support
target_link_libraries(ht talyfem)

# automatically generate documentation with Doxygen when running 'make docs'
add_doxygen(TARGET ht INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/include)

# generate reproducibility information header (build_info_ssht.h)
add_taly_reproducibility(TARGET ht)

# add tests
enable_testing()
foreach(test_idx RANGE 0 53)
  add_test(NAME "tutorial_ht_${test_idx}"
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_ht.py ${test_idx}
    )
  set_property(TEST "tutorial_ht_${test_idx}"
    PROPERTY ENVIRONMENT
    EXECUTABLE_PATH=$<TARGET_FILE:ht>
    TEST_DATA_DIR=${CMAKE_CURRENT_SOURCE_DIR}/tests
    TALYFEM_DIR=$<TARGET_PROPERTY:talyfem,SOURCE_DIR>
    )
endforeach()
