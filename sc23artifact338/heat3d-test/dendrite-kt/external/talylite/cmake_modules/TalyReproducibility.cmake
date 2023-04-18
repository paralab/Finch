include (CMakeParseArguments)

find_package(PythonInterp 2.7 REQUIRED)

function(add_taly_reproducibility)
  set(oneValueArgs TARGET BUILD_INFO_GEN_PATH)
  set(multiValueArgs DEPENDS)
  cmake_parse_arguments(REPRO "" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  if (NOT REPRO_TARGET)
    message(FATAL_ERROR "add_taly_reproducibility: must specify TARGET")
  endif()

  if (NOT REPRO_DEPENDS)
    get_target_property(REPRO_DEPENDS ${REPRO_TARGET} SOURCES)
  endif()

  if (NOT REPRO_BUILD_INFO_GEN_PATH)
    if(EXISTS "${CMAKE_SOURCE_DIR}/python_scripts/build_info_gen.py")
      set(REPRO_BUILD_INFO_GEN_PATH "${CMAKE_SOURCE_DIR}/python_scripts/build_info_gen.py")
    elseif(EXISTS "${TALYFEM_PYTHON_SCRIPTS}/build_info_gen.py")
      set(REPRO_BUILD_INFO_GEN_PATH "${TALYFEM_PYTHON_SCRIPTS}/build_info_gen.py")
    endif()
  endif()

  if (NOT EXISTS "${REPRO_BUILD_INFO_GEN_PATH}")
    message(FATAL_ERROR "add_taly_reproducibility: missing build_info_gen.py (${REPRO_BUILD_INFO_GEN_PATH})")
  endif()

  set(REPRO_OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/build_info_${REPRO_TARGET}.h)

  add_custom_command(
    OUTPUT ${REPRO_OUTPUT}
    COMMAND ${CMAKE_COMMAND} -E env
        CMAKE_CACHEFILE=${CMAKE_BINARY_DIR}/CMakeCache.txt
        CMAKE_CXX_COMPILER_ID="${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}"
        ${PYTHON_EXECUTABLE} ${REPRO_BUILD_INFO_GEN_PATH} --output ${REPRO_OUTPUT} --struct-name BuildInfo_${REPRO_TARGET}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Generating reproducibility information for ${REPRO_TARGET}..."
    DEPENDS ${CMAKE_BINARY_DIR}/CMakeCache.txt ${REPRO_DEPENDS}
    VERBATIM
  )

  # Add the reproducibility file to the build
  target_sources(${REPRO_TARGET} PRIVATE ${REPRO_OUTPUT})
  target_include_directories(${REPRO_TARGET} PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
endfunction()
