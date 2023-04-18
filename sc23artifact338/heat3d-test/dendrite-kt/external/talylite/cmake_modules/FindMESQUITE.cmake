find_path(MESQUITE_INCLUDE_DIR
  NAMES Mesquite.hpp
  HINTS ${MESQUITE_DIR} $ENV{MESQUITE_DIR}
  PATH_SUFFIXES include
)

find_library(MESQUITE_LIBRARY
  NAMES mesquite
  HINTS ${MESQUITE_DIR} $ENV{MESQUITE_DIR}
  PATH_SUFFIXES lib
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MESQUITE
  REQUIRED_VARS MESQUITE_LIBRARY MESQUITE_INCLUDE_DIR
)
