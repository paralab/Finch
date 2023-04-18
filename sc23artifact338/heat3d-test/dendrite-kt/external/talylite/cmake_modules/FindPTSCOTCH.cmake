# look for ptscotch.h
find_path(PTSCOTCH_INCLUDE_DIRS
  NAMES ptscotch.h
  HINTS ${SCOTCH_DIR}
  PATH_SUFFIXES include
)

# look for the compiled scotch and ptscotch library files
find_library(PTSCOTCH_LIBRARY_SCOTCH
  NAMES scotch
  HINTS ${SCOTCH_DIR}
  PATH_SUFFIXES lib
)
find_library(PTSCOTCH_LIBRARY_PTSCOTCH
  NAMES ptscotch
  HINTS ${SCOTCH_DIR}
  PATH_SUFFIXES lib
)

# PTSCOTCH offers two error handlers, one which prints errors and continues,
# and one which prints errors and aborts. Switching between the libraries is
# done by linking between one of two libraries that contain the error handler
# functions.
# By default, this script will use the first way (report and continue).
# Set PTSCOTCH_ERREXIT to true to use the second way (abort on error).
if(PTSCOTCH_ERREXIT)
  find_library(PTSCOTCH_LIBRARY_PTSCOTCHERR
    NAMES ptscotcherrexit
    HINTS ${SCOTCH_DIR}
    PATH_SUFFIXES lib
  )
else()
  find_library(PTSCOTCH_LIBRARY_PTSCOTCHERR
    NAMES ptscotcherr
    HINTS ${SCOTCH_DIR}
    PATH_SUFFIXES lib
  )
endif()

# order is important: scotch must get linked after ptscotch
set(PTSCOTCH_LIBRARIES
  ${PTSCOTCH_LIBRARY_PTSCOTCH}
  ${PTSCOTCH_LIBRARY_SCOTCH}
  ${PTSCOTCH_LIBRARY_PTSCOTCHERR}
)

# PTSCOTCH optionally includes a ParMETIS shim for a subset of ParMETIS functions.
find_library(PTSCOTCH_LIBRARY_PTSCOTCHPARMETIS
  NAMES ptscotchparmetis
  HINTS ${SCOTCH_DIR}
  PATH_SUFFIXES lib
)

# only link to it if the "ParMETIS" component was requested
if (PTSCOTCH_FIND_COMPONENTS AND ";${PTSCOTCH_FIND_COMPONENTS};" MATCHES ";ParMETIS;")
  # if the component was found, add it to the libraries to link against
  # order is important: must go *before* the main scotch libraries
  if(PTSCOTCH_LIBRARY_PTSCOTCHPARMETIS)
    set(PTSCOTCH_LIBRARIES ${PTSCOTCH_LIBRARY_PTSCOTCHPARMETIS} ${PTSCOTCH_LIBRARIES} )
  else()
    # otherwise, raise an error
    message(SEND_ERROR "PTSCOTCH ParMETIS shim is missing (libptscotchparmetis).")
  endif()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PTSCOTCH
  REQUIRED_VARS PTSCOTCH_LIBRARIES PTSCOTCH_INCLUDE_DIRS
)