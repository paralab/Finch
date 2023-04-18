include (CMakeParseArguments)

find_package(Doxygen)

if(NOT TARGET docs)
  add_custom_target(docs)
endif()

function (add_doxygen)
  set(noValueArgs GENERATE_TAGFILE)
  set(oneValueArgs TARGET DOXYGEN_BASE)
  set(multiValueArgs INCLUDE_DIRS EXCLUDE)
  cmake_parse_arguments(add_doxygen "${noValueArgs}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  if (NOT add_doxygen_TARGET)
    message(FATAL_ERROR "add_doxygen: must specify TARGET")
  endif()

  if (NOT add_doxygen_INCLUDE_DIRS)
    message(SEND_ERROR "add_doxygen: must specify INCLUDE_DIRS")
  endif()

  # Variables used in doxygen.cfg.in (convert from semicolon-separated to space-separated)
  string(REPLACE ";" " " DOXYGEN_INPUT "${add_doxygen_INCLUDE_DIRS}")
  string(REPLACE ";" " " DOXYGEN_EXCLUDE "${add_doxygen_EXCLUDE}")
  set(GENERATE_TAGFILE "docs/${add_doxygen_TARGET}.tag")

  if (NOT add_doxygen_DOXYGEN_BASE)
    if(EXISTS "${CMAKE_SOURCE_DIR}/docs/doxygen.cfg.in")
      set(add_doxygen_DOXYGEN_BASE "${CMAKE_SOURCE_DIR}/docs/doxygen.cfg.in")
    elseif(EXISTS "${TALYFEM_DOXYGEN_BASE}")
      set(add_doxygen_DOXYGEN_BASE "${TALYFEM_DOXYGEN_BASE}")
    endif()
  endif()

  if (NOT EXISTS "${add_doxygen_DOXYGEN_BASE}")
    message(FATAL_ERROR "add_doxygen: DOXYGEN_BASE not found (${add_doxygen_DOXYGEN_BASE})")
  endif()

  if (DOXYGEN_FOUND)
    configure_file(${add_doxygen_DOXYGEN_BASE}
      ${CMAKE_CURRENT_BINARY_DIR}/doxygen.cfg @ONLY)

    add_custom_target(${add_doxygen_TARGET}-docs
      COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/doxygen.cfg
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      COMMENT "Generating documentation with Doxygen for ${add_doxygen_TARGET}..."
      DEPENDS ${add_doxygen_FILES}
    )
  else()
    add_custom_target(${add_doxygen_TARGET}-docs
      COMMAND ${CMAKE_COMMAND} -E echo "Doxygen is not installed. Cannot generate docs."
    )
  endif()

  add_dependencies(docs ${add_doxygen_TARGET}-docs)
endfunction()