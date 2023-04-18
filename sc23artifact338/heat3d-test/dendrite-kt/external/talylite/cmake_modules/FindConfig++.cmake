find_path(CONFIG++_INCLUDE_DIR
	NAMES libconfig.h++
	HINTS ${LIBCONFIG_DIR} $ENV{LIBCONFIG_DIR}
	PATH_SUFFIXES include
)

find_library(CONFIG++_LIBRARY
	NAMES config++
	HINTS ${LIBCONFIG_DIR} $ENV{LIBCONFIG_DIR}
	PATH_SUFFIXES lib
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Config++
  REQUIRED_VARS CONFIG++_LIBRARY CONFIG++_INCLUDE_DIR
)
