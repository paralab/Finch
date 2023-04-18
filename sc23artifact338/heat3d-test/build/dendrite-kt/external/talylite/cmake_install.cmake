# Install script for directory: /work/08517/heisler/frontera/dendritenew/dendrite-kt/external/talylite

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/opt/apps/gcc/8.3.0/bin/objdump")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/work/08517/heisler/frontera/dendritenew/build/dendrite-kt/external/talylite/libtalyfem.a")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE DIRECTORY FILES "/work/08517/heisler/frontera/dendritenew/dendrite-kt/external/talylite/include/talyfem")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nanoflann" TYPE FILE FILES "/work/08517/heisler/frontera/dendritenew/dendrite-kt/external/talylite/external/nanoflann/include/nanoflann.hpp")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/talyfem/cmake" TYPE FILE FILES
    "/work/08517/heisler/frontera/dendritenew/dendrite-kt/external/talylite/cmake_modules/TalyDoxygen.cmake"
    "/work/08517/heisler/frontera/dendritenew/dendrite-kt/external/talylite/cmake_modules/TalyReproducibility.cmake"
    "/work/08517/heisler/frontera/dendritenew/dendrite-kt/external/talylite/cmake_modules/TalyStaticAnalysis.cmake"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/talyfem/docs" TYPE FILE FILES "/work/08517/heisler/frontera/dendritenew/dendrite-kt/external/talylite/docs/doxygen.cfg.in")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/talyfem/python_scripts/" TYPE DIRECTORY FILES "/work/08517/heisler/frontera/dendritenew/dendrite-kt/external/talylite/python_scripts/" REGEX "/[^/]*\\.pyc$" EXCLUDE)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/share/talyfem/cmake/talyfemConfig.cmake")
    file(DIFFERENT _cmake_export_file_changed FILES
         "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/share/talyfem/cmake/talyfemConfig.cmake"
         "/work/08517/heisler/frontera/dendritenew/build/dendrite-kt/external/talylite/CMakeFiles/Export/e5bccdd3e74ece9bd4514f17a732d09b/talyfemConfig.cmake")
    if(_cmake_export_file_changed)
      file(GLOB _cmake_old_config_files "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/share/talyfem/cmake/talyfemConfig-*.cmake")
      if(_cmake_old_config_files)
        string(REPLACE ";" ", " _cmake_old_config_files_text "${_cmake_old_config_files}")
        message(STATUS "Old export file \"$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/share/talyfem/cmake/talyfemConfig.cmake\" will be replaced.  Removing files [${_cmake_old_config_files_text}].")
        unset(_cmake_old_config_files_text)
        file(REMOVE ${_cmake_old_config_files})
      endif()
      unset(_cmake_old_config_files)
    endif()
    unset(_cmake_export_file_changed)
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/talyfem/cmake" TYPE FILE FILES "/work/08517/heisler/frontera/dendritenew/build/dendrite-kt/external/talylite/CMakeFiles/Export/e5bccdd3e74ece9bd4514f17a732d09b/talyfemConfig.cmake")
  if(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Rr][Ee][Ll][Ee][Aa][Ss][Ee])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/talyfem/cmake" TYPE FILE FILES "/work/08517/heisler/frontera/dendritenew/build/dendrite-kt/external/talylite/CMakeFiles/Export/e5bccdd3e74ece9bd4514f17a732d09b/talyfemConfig-release.cmake")
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(APPEND /usr/local/share/talyfem/cmake/talyfemConfig.cmake
    "list(APPEND CMAKE_MODULE_PATH /usr/local/share/talyfem/cmake/)
    set(TALYFEM_DOXYGEN_BASE /usr/local/share/talyfem/docs/doxygen.cfg.in)
    set(TALYFEM_PYTHON_SCRIPTS /usr/local/share/talyfem/python_scripts/)")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/work/08517/heisler/frontera/dendritenew/build/dendrite-kt/external/talylite/external/b64/cmake_install.cmake")

endif()

