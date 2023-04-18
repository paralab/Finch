include (CMakeParseArguments)

# Try and find clang-tidy
find_program(CLANG_TIDY_EXECUTABLE
  NAMES clang-tidy clang-tidy-3.8 clang-tidy-3.9 clang-tidy-4.0
  DOC "Path to clang-tidy executable"
)

if (CLANG_TIDY_EXECUTABLE)
  if (CMAKE_VERSION VERSION_LESS 3.7)
    set(CLANG_TIDY_EXECUTABLE "")
    message(STATUS "CMake version less than 3.7, which does not support clang-tidy.")
  else()
    message(STATUS "Found clang-tidy: ${CLANG_TIDY_EXECUTABLE}")
  endif()
else()
  message(STATUS "clang-tidy not found.")
endif()

function (add_taly_static_analysis)
  set(noValueArgs WARNINGS_AS_ERRORS)  # TODO
  set(oneValueArgs TARGET)
  set(multiValueArgs CHECKS)
  cmake_parse_arguments(STATIC_ANALYSIS "" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  if (NOT STATIC_ANALYSIS_TARGET)
    message(FATAL_ERROR "add_taly_static_analysis: must specify TARGET")
  endif()

  if (NOT STATIC_ANALYSIS_CHECKS)
    set(STATIC_ANALYSIS_CHECKS
      -readability-simplify-boolean-expr
      -clang-analyzer-unix.API  # skip checks for malloc(0)
      google-*
      -google-readability-casting
      -google-readability-braces-around-statements
      -google-readability-todo)
  endif()

  string(REPLACE ";" "," STATIC_ANALYSIS_CHECKS "${STATIC_ANALYSIS_CHECKS}")
  set(DO_CLANG_TIDY "${CLANG_TIDY_EXECUTABLE}" "-checks=${STATIC_ANALYSIS_CHECKS}")

  if (CLANG_TIDY_EXECUTABLE)
    set_target_properties(${STATIC_ANALYSIS_TARGET} PROPERTIES
      CXX_CLANG_TIDY "${DO_CLANG_TIDY}"
    )
  endif()
endfunction()