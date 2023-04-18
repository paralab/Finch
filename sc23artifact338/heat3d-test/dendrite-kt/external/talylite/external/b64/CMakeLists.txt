cmake_minimum_required(VERSION 3.5.0 FATAL_ERROR)

set(B64_INCLUDES
  include/cdecode.h
  include/cencode.h
  include/decode.h
  include/encode.h
)

set(B64_SOURCES
  src/cdecode.cpp
  src/cencode.cpp
)

# Compile as "OBJECT" library so it gets included inside libtalyfem.a
# From here: https://cmake.org/Wiki/CMake/Tutorials/Object_Library
add_library(b64 OBJECT
  ${B64_INCLUDES}
  ${B64_SOURCES}
)

target_include_directories(b64
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include
)
