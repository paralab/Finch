cmake_minimum_required(VERSION 3.5.0 FATAL_ERROR)

cmake_host_system_information(RESULT HOSTNAME QUERY HOSTNAME)

# If no compiler is specified, try using the `mpiicpc` compiler when running
# on CyEnce or Condo.
# This must be done before project() is called.
if (NOT DEFINED CMAKE_CXX_COMPILER AND ((HOSTNAME MATCHES "^condo\\d*") OR (HOSTNAME MATCHES "^cyence.*")))
  find_program(FOUND_MPIICPC "mpiicpc")
  if (FOUND_MPIICPC)
    message("[CyEnce/Condo Intel compiler detected] Overriding CXX: ${FOUND_MPIICPC}")
    set(CMAKE_CXX_COMPILER "${FOUND_MPIICPC}")
  endif()
endif()

if (NOT DEFINED CMAKE_CXX_COMPILER AND HOSTNAME MATCHES "^comet")
  find_program(FOUND_MPICXX "mpicxx")
  if (FOUND_MPICXX)
    message("[Comet Intel compiler detected] Overriding CXX: ${FOUND_MPICXX}")
    set(CMAKE_CXX_COMPILER "${FOUND_MPICXX}")
  endif()
endif()

project(talyfem)

# Add the cmake_modules directory to the CMake module search path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules")

# Try to find PETSc 3.6+
find_package(PETSc 3.6 COMPONENTS C REQUIRED)

# Try to find MPI
if (NOT DEFINED MPI_HOME)
  set(MPI_HOME ${PETSC_DIR}/${PETSC_ARCH})
endif()
find_package(MPI COMPONENTS C REQUIRED)
string(STRIP "${MPI_C_LINK_FLAGS}" MPI_C_LINK_FLAGS)

# Try to find ParMETIS
find_package(ParMETIS REQUIRED)

# Try to find libconfig (C++ version)
find_package(Config++)
if (NOT ${CONFIG++_FOUND})
  # TODO - build bundled libconfig
  message(FATAL_ERROR "Libconfig++ could not be located.")
endif()

# Add external/b64/CMakeLists.txt (and the b64 object target)
add_subdirectory(external/b64)

# Define all the include and source files for TalyFEM
set(TALYFEM_HEADERS
  include/talyfem/integrator/function_integrator.h
  include/talyfem/integrator/value_function.h

  include/talyfem/grid/segment.h
  include/talyfem/grid/grid_types/gridbox2d.h
  include/talyfem/grid/grid_types/grid1d.h
  include/talyfem/grid/grid_types/grid.h
  include/talyfem/grid/grid_types/grid2d.h
  include/talyfem/grid/grid_types/gridbox3d.h
  include/talyfem/grid/cremoteneighbor.h
  include/talyfem/grid/surfaceindicator.h
  include/talyfem/grid/kdtree.h
  include/talyfem/grid/grid_common.h
  include/talyfem/grid/node.h
  include/talyfem/grid/elem.h
  include/talyfem/grid/elem-types.h
  include/talyfem/grid/nodedata.h
  include/talyfem/grid/elem_common.h
  include/talyfem/grid/gridfield.h
  include/talyfem/grid/femelm.h
  include/talyfem/grid/zeroptv.h
  include/talyfem/grid/sharedetails.h
  include/talyfem/grid/grid-types.h
  include/talyfem/grid/shareinfo.h
  include/talyfem/grid/nodeindicator.h
  include/talyfem/grid/elem_types/elem3dhexahedral.h
  include/talyfem/grid/elem_types/elem2dbox.h
  include/talyfem/grid/elem_types/elem3dtetrahedral.h
  include/talyfem/grid/elem_types/elem2dtriangle.h
  include/talyfem/grid/mesquite.h
  include/talyfem/grid/nodeid_types.h
        include/talyfem/grid/elem-quality.h

  include/talyfem/basis/constants.h
  include/talyfem/basis/elemnodes.h
  include/talyfem/basis/box_linear/box_linear_impl.h
  include/talyfem/basis/box_linear/box_linear.h
  include/talyfem/basis/basis.h
  include/talyfem/basis/box_cubic/box_cubic_impl.h
  include/talyfem/basis/box_cubic/box_cubic.h
  include/talyfem/basis/nonlinear.h
  include/talyfem/basis/point.h
  include/talyfem/basis/tri_linear/tri_linear_impl.h
  include/talyfem/basis/tri_linear/tri_linear.h
  include/talyfem/basis/tri_quadratic/tri_quadratic_impl.h
  include/talyfem/basis/tri_quadratic/tri_quadratic.h
  include/talyfem/basis/tri_cubic_10pts/tri_cubic10_impl.h
  include/talyfem/basis/tri_cubic_10pts/tri_cubic10.h
  include/talyfem/basis/tet_linear/tet_linear_impl.h
  include/talyfem/basis/tet_linear/tet_linear.h
  include/talyfem/basis/tet_quadratic/tet_quadratic_impl.h
  include/talyfem/basis/tet_quadratic/tet_quadratic.h
  include/talyfem/basis/tet_cubic/tet_cubic_impl.h
  include/talyfem/basis/tet_cubic/tet_cubic.h
  include/talyfem/basis/box_hermite/box_hermite.h
  include/talyfem/basis/box_hermite/box_hermite_impl.h
  include/talyfem/basis/matrix/matrix.h
  include/talyfem/basis/matrix/matrix_ops.h
  include/talyfem/basis/itg_pts/itg_pts.h
  include/talyfem/basis/itg_pts/box_2d.h
  include/talyfem/basis/itg_pts/box_1d.h
  include/talyfem/basis/itg_pts/box_3d.h
  include/talyfem/basis/itg_pts/tet.h
  include/talyfem/basis/itg_pts/box_surfaces.h
  include/talyfem/basis/itg_pts/tri_2d.h
  include/talyfem/basis/basis_common.h
  include/talyfem/basis/type_traits.h
  include/talyfem/basis/box_quadratic/box_quadratic_impl.h
  include/talyfem/basis/box_quadratic/box_quadratic.h
  include/talyfem/basis/mat3.h

  include/talyfem/fem/periodic_exchanger.h
  include/talyfem/fem/boundary_conditions.h
  include/talyfem/fem/essential_boundary_condition.h
  include/talyfem/fem/periodic_bounds.h
  include/talyfem/fem/preallocator.h
  include/talyfem/fem/periodic_data.h
  include/talyfem/fem/zero_boundary_condition.h
  include/talyfem/fem/preallocator_perfect.h
  include/talyfem/fem/base_boundary_condition.h
  include/talyfem/fem/cequation.h
  include/talyfem/fem/preallocator_original.h

  include/talyfem/stabilizer/tezduyar_upwind.h

  include/talyfem/domain_decomposition/mesh_partition.h
  include/talyfem/domain_decomposition/scotch_wrapper.h
  include/talyfem/domain_decomposition/cross_domain_communicator.h

  include/talyfem/utils/utils.h
  include/talyfem/utils/timers.h
  include/talyfem/utils/print_utils_variadic.h
  include/talyfem/utils/test_utils.h
  include/talyfem/utils/reproducibility.h
  include/talyfem/utils/macros.h

  include/talyfem/file_io/tecplot_dd.h
  include/talyfem/file_io/surface_io.h
  include/talyfem/file_io/tecplot_grid.h
  include/talyfem/file_io/file_io.h
  include/talyfem/file_io/h5_io.h
  include/talyfem/file_io/gmsh_grid.h
  include/talyfem/file_io/tecplot_gf.h
  include/talyfem/file_io/tecplot_ascii.h
  include/talyfem/file_io/h5_grid.h
  include/talyfem/file_io/tecplot_io.h
  include/talyfem/file_io/h5_gf.h
  include/talyfem/file_io/common.h
  include/talyfem/file_io/gmsh_io.h

  include/talyfem/common/indexer2d.h
  include/talyfem/common/comm.h
  include/talyfem/common/exceptions.h
  include/talyfem/common/pack_comm.h
  include/talyfem/common/petsc_logging.h
  include/talyfem/common/indexer3d.h

  include/talyfem/data_structures/constexpr_array.h
  include/talyfem/data_structures/zeromatrix.h
  include/talyfem/data_structures/matrix4.h
  include/talyfem/data_structures/zeroarray.h

  include/talyfem/input_data/input_data.h
        include/talyfem/utils/expressions.h

  )

set(TALYFEM_SOURCES
  src/basis/basis.cpp
  src/basis/itg_pts/itg_pts.cpp

  src/common/comm.cpp
  src/common/exceptions.cpp
  src/common/pack_comm.cpp
  src/common/petsc_logging.cpp

  src/domain_decomposition/cross_domain_communicator.cpp
  src/domain_decomposition/mesh_partition.cpp

  src/fem/base_boundary_condition.cpp
  src/fem/essential_boundary_condition.cpp
  src/fem/periodic_bounds.cpp
  src/fem/periodic_data.cpp
  src/fem/periodic_exchanger.cpp
  src/fem/preallocator_original.cpp
  src/fem/zero_boundary_condition.cpp

  src/file_io/gmsh_grid.cpp
  src/file_io/gmsh_io.cpp
  src/file_io/h5_grid.cpp
  src/file_io/h5_io.cpp
  src/file_io/surface_io.cpp
  src/file_io/tecplot_ascii.cpp
  src/file_io/tecplot_grid.cpp

  src/grid/cremoteneighbor.cpp
  src/grid/elem.cpp
  src/grid/elem_types.cpp
  src/grid/elem_types/elem1d.cpp
  src/grid/elem_types/elem2dbox.cpp
  src/grid/elem_types/elem2dtriangle.cpp
  src/grid/elem_types/elem3dhexahedral.cpp
  src/grid/elem_types/elem3dtetrahedral.cpp
  src/grid/femelm.cpp
  src/grid/grid_types/grid1d.cpp
  src/grid/grid_types/grid2d.cpp
  src/grid/grid_types/gridbox2d.cpp
  src/grid/grid_types/gridbox3d.cpp
  src/grid/grid_types/grid.cpp
  src/grid/kdtree.cpp
  src/grid/mesquite.cpp
  src/grid/node.cpp
  src/grid/segment.cpp
  src/grid/shareinfo.cpp
  src/grid/surfaceindicator.cpp
  src/grid/zeroptv.cpp

  src/input_data/input_data.cpp

  src/integrator/function_integrator.cpp

  src/math/math.cpp

  src/stabilizer/tezduyar_upwind.cpp

        src/utils/expressions.cpp
  src/utils/reproducibility.cpp
  src/utils/test_utils.cpp
  src/utils/timers.cpp
  src/utils/utils.cpp
)
option(ENABLE_4D "Enable 4D Support" OFF)

if (ENABLE_4D)
    list(APPEND TALYFEM_SOURCES
            "src/grid/grid_types/gridbox4d.cpp"
            "src/grid/elem_types/elem4dtesseract.cpp"
            "src/basis/itg_pts/box_surfaces.cpp")
    list(APPEND TALYFEM_HEADERS
            "include/talyfem/grid/grid_types/gridbox4d.h")
endif ()
if(ENABLE_4D)
    message("4D computation ON")
    endif()
# Create the talyfem library
add_library(talyfem STATIC
  ${TALYFEM_SOURCES}
  ${TALYFEM_HEADERS}

  # libb64 is included as an object library, so we must include it here instead of in target_link_libraries()
  $<TARGET_OBJECTS:b64>
)

# Require the compiler to enable C++11 support (or above...)
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel" AND CMAKE_VERSION VERSION_LESS 3.6)
  # CMake didn't get support for target_compile_features with the Intel compiler until CMake 3.6

  # Intel compiler versions 14 and below have problems compiling (due to incomplete c++11 support)
  if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 15)
    message(SEND_ERROR "TalyFEM requires Intel compiler version 15 or above.")
  endif()

  # so we just manually set the c++11 compiler flags when using older versions of CMake
  target_compile_options(talyfem PUBLIC "-std=c++11")
else()
  # On other compilers (and Intel compilers with CMake 3.6+), we can use target_compile_features
  # (this is the recommended way of activating c++11 support)
  target_compile_features(talyfem PUBLIC cxx_override cxx_auto_type cxx_variadic_templates cxx_constexpr)
endif()


# On Linux machines we need to link against librt for clock_gettime
# (well, *sometimes*, depending on the verison of libc), which is  used
# in the HiRestimer class. If we are using an up-to-date version of libc,
# it won't hurt to link against librt (we just won't consume any symbols).
# OSX does not have this library, so we must not include it on APPLE systems.
if (UNIX AND NOT APPLE)
  target_link_libraries(talyfem PUBLIC rt)
endif()

# List include directories talyfem (and users) will need to include
target_include_directories(talyfem
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include>
  PUBLIC ${PETSC_INCLUDES}
  PUBLIC ${MPI_C_INCLUDE_PATH}
  PUBLIC ${CONFIG++_INCLUDE_DIR}
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/external/nanoflann/include> $<INSTALL_INTERFACE:include/nanoflann>
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/external/b64/include
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/external/exprtk/include

)

# List other libraries talyfem users will need to link against
target_link_libraries(talyfem
  LINK_PUBLIC ${PETSC_LIBRARIES} ${MPI_C_LIBRARIES} ${MPI_C_LINK_FLAGS} ${PARMETIS_LIBRARIES} ${CONFIG++_LIBRARY}
)

option(ENABLE_HDF5 "Enable HDF5 file I/O (requires HDF5 library)" OFF)
if (ENABLE_HDF5)
  set(HDF5_PREFER_PARALLEL TRUE)
  find_package(HDF5 REQUIRED COMPONENTS C)
  if (NOT HDF5_IS_PARALLEL)
    message(WARN "HDF5 library is not parallel. You may get run-time errors using HDF5.")
  endif()

  # add HDF5 to the TalyFEM build configuration
  target_compile_definitions(talyfem PUBLIC ${HDF5_DEFINITIONS} -DENABLE_HDF5)
  target_include_directories(talyfem PUBLIC ${HDF5_INCLUDE_DIRS})
  target_link_libraries(talyfem LINK_PUBLIC ${HDF5_LIBRARIES})
endif()


if (ENABLE_4D)
    target_compile_definitions(talyfem PUBLIC ${HDF5_DEFINITIONS} -DENABLE_4D)
endif ()

option(ENABLE_SCOTCH "Use PTSCOTCH for domain decomposition partition instead of ParMETIS (requires the PTSCOTCH library) (still requires ParMETIS to build the dual graph)" OFF)
if (ENABLE_SCOTCH)
  find_package(PTSCOTCH REQUIRED COMPONENTS ParMETIS)

  # add PTSCOTCH to the TalyFEM build configuration
  target_compile_definitions(talyfem PRIVATE -DENABLE_SCOTCH)
  target_include_directories(talyfem PRIVATE ${PTSCOTCH_INCLUDE_DIRS})
  target_link_libraries(talyfem LINK_PUBLIC ${PTSCOTCH_LIBRARIES})
endif()


option(ENABLE_MESQUITE "Compile TalyFEM-Mesquite mesh adapter (requires the Mesquite library)" OFF)
if (ENABLE_MESQUITE)
  find_package(MESQUITE REQUIRED)

  # add MESQUITE to the TalyFEM build configuration
  target_compile_definitions(talyfem PUBLIC -DENABLE_MESQUITE)
  target_include_directories(talyfem PUBLIC ${MESQUITE_INCLUDE_DIR})
  target_link_libraries(talyfem LINK_PUBLIC ${MESQUITE_LIBRARY})

  # add the mesquite test
  add_subdirectory(tests/mesquite)
endif()

# The gauss point code uses constexpr std::array access. Some systems' STL
# implementations do not support this (even though the compiler can handle it).
# There is a replacement for std::array  for this case, but it must be
# explicitly turned on. We try to compile a file that uses constexpr std::array
# access to see if we need to enable our override.

include (CMakePushCheckState)
include (CheckCXXSourceCompiles)

cmake_push_check_state()
set(CMAKE_REQUIRED_FLAGS ${CMAKE_REQUIRED_FLAGS} ${CMAKE_CXX11_STANDARD_COMPILE_OPTION})
check_cxx_source_compiles("
#include <array>
constexpr std::array<double, 3> arr = {{ 1.0, 2.0, 3.0 }};
constexpr double test = arr[1];
int main() {
  return (test != 2.0);
}
" SUPPORTS_CONSTEXPR_STDARRAY)
cmake_pop_check_state()

if (NOT SUPPORTS_CONSTEXPR_STDARRAY)
  message("Compiler does not appear to support constexpr std::array accesses - using replacement.")
  target_compile_definitions(talyfem PUBLIC -DNO_STDARRAY)
endif()

# Add reproducibility information
include (TalyReproducibility)
set(TALYFEM_PYTHON_SCRIPTS "${CMAKE_CURRENT_SOURCE_DIR}/python_scripts")
add_taly_reproducibility(TARGET talyfem DEPENDS ${TALYFEM_SOURCES} ${TALYFEM_HEADERS})

# Add automatic doxygen generation
include (TalyDoxygen)
set(TALYFEM_DOXYGEN_BASE "${CMAKE_CURRENT_SOURCE_DIR}/docs/doxygen.cfg.in")
add_doxygen(TARGET talyfem INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Add automatic static analysis checks with clang-tidy during compilation
include (TalyStaticAnalysis)
add_taly_static_analysis(TARGET talyfem)

# Make 'make install' install the correct files
install(TARGETS talyfem
  EXPORT talyfemConfig
  ARCHIVE DESTINATION lib  # relative dir to put .a files in
  LIBRARY DESTINATION lib  # relative dir to put .so files in
  RUNTIME DESTINATION bin  # relative dir to put executables in
)

# also copy the entire include/ directory and any external dependencies we depend on
install(DIRECTORY include/talyfem DESTINATION include)
install(FILES external/nanoflann/include/nanoflann.hpp DESTINATION include/nanoflann/)

# also install some helpful CMake modules for other projects to use
install(FILES
    cmake_modules/TalyDoxygen.cmake
    cmake_modules/TalyReproducibility.cmake
    cmake_modules/TalyStaticAnalysis.cmake
    DESTINATION share/talyfem/cmake/)

# also install the Doxygen config file for sub projects using TalyDoxygen
install(FILES
    docs/doxygen.cfg.in
    DESTINATION share/talyfem/docs/)

# and the build info generator for reproducibility info for sub projects
install(DIRECTORY python_scripts/ DESTINATION share/talyfem/python_scripts/
    PATTERN "*.pyc" EXCLUDE)

# finally, export the TalyFEM target to share/talyfem/talyfemConfig.cmake
install(EXPORT talyfemConfig DESTINATION share/talyfem/cmake/)

# this runs a piece of CMake code after exporting talyfemConfig.cmake during install
# that adds share/talyfem/talyfemConfig.cmake's directory to the module search path,
# so that users can directly include() the files from cmake_modules we installed above
# this only happens during `make install`
# if you know a better way to make the taly modules automatically discoverable
# after find_package(), please improve this
install(CODE "file(APPEND ${CMAKE_INSTALL_PREFIX}/share/talyfem/cmake/talyfemConfig.cmake
    \"list(APPEND CMAKE_MODULE_PATH ${CMAKE_INSTALL_PREFIX}/share/talyfem/cmake/)
    set(TALYFEM_DOXYGEN_BASE ${CMAKE_INSTALL_PREFIX}/share/talyfem/docs/doxygen.cfg.in)
    set(TALYFEM_PYTHON_SCRIPTS ${CMAKE_INSTALL_PREFIX}/share/talyfem/python_scripts/)\")")

# Also manually export the talyfem target to a talyfemConfig.cmake in the build dir
# that can be  imported by other CMake projects directly from the build directory
# (no install needed)
# we do similar kung-fu as for install() above: we manually create talyfemConfig.cmake
# first and write a command to append to CMAKE_MODULE_PATH so users can include()
# the files in cmake_modules/, *then* append the exported target to that file.
# if you know a better way to make the taly modules automatically discoverable
# after find_package(), please improve this
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/talyfemConfig.cmake
    "list(APPEND CMAKE_MODULE_PATH \"${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules\")
    set(TALYFEM_DOXYGEN_BASE \"${CMAKE_CURRENT_SOURCE_DIR}/docs/doxygen.cfg.in\")
    set(TALYFEM_PYTHON_SCRIPTS \"${CMAKE_CURRENT_SOURCE_DIR}/python_scripts\")\n\n")
export(TARGETS talyfem APPEND FILE talyfemConfig.cmake)

# We must call enable_testing() here, not just in the subdirectories, in order
# for tests to show up in the top-level `make test` target. If this was left
# out, we would have to cd into each sub-directory and run `make test` on
# each individual tutorial/test.
enable_testing()

## Add the tutorial and test directories to the build
#option(TALYFEM_BUILD_TUTORIALS "Build tutorial applications (HT, BT)" ON)
#if (TALYFEM_BUILD_TUTORIALS)
#    add_subdirectory(tutorials/steadystate_heat)
#    add_subdirectory(tutorials/transient_heat)
#    add_subdirectory(tutorials/bratu)
#endif()
#
#option(TALYFEM_BUILD_TESTS "Build unit tests" ON)
#if (TALYFEM_BUILD_TESTS)
#    add_subdirectory(tests/basic_io)
#    add_subdirectory(tests/basis)
#    add_subdirectory(tests/cross_domain_comm)
#    add_subdirectory(tests/element)
#    add_subdirectory(tests/integrator)
#    add_subdirectory(tests/periodic_dd_map)
#    add_subdirectory(tests/periodic_fill)
#    add_subdirectory(tests/stabilizer)
#    add_subdirectory(tests/tecplot_io)
#endif()
