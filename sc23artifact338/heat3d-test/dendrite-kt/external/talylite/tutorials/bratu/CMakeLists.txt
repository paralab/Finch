cmake_minimum_required(VERSION 3.5.0 FATAL_ERROR)

project(bratu)

add_executable(bt
  # header files (not needed, but added so they show up in IDEs)
  include/BTEquation.h
  include/BTGridField.h
  include/BTInputData.h
  include/BTNodeData.h

  # source files to compile
  src/main.cpp
)

# add the ./include directory to the include search path so we can #include
# the header files we put in include/
target_include_directories(bt PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

if (NOT TARGET talyfem)
  find_package(talyfem REQUIRED)
  include(TalyDoxygen)  # for add_doxygen()
  include(TalyReproducibility)  # for add_taly_reproducibility()
endif()

# this will also copy all the compiler definitions, etc. that TalyFEM was
# compiled with (that were marked as PUBLIC) - e.g. this will automatically
# add "-DENABLE_HDF5" if TalyFEM was compiled with HDF5 support
target_link_libraries(bt talyfem)

# automatically generate documentation with Doxygen when running 'make docs'
add_doxygen(TARGET bt INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/include)

# generate reproducibility information header (build_info_ssht.h)
add_taly_reproducibility(TARGET bt)

enable_testing()
foreach(test_idx RANGE 0 29)
  add_test(NAME "tutorial_bt_${test_idx}"
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_bratu.py ${test_idx}
    )
  set_property(TEST "tutorial_bt_${test_idx}"
    PROPERTY ENVIRONMENT
    EXECUTABLE_PATH=$<TARGET_FILE:bt>
    TEST_DATA_DIR=${CMAKE_CURRENT_SOURCE_DIR}/tests
    TALYFEM_DIR=$<TARGET_PROPERTY:talyfem,SOURCE_DIR>
    )
endforeach()
