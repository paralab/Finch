## CMAKE for Dendro-KT
cmake_minimum_required(VERSION 2.8)
project(Dendro-KT)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)

# For now we just make it compulsory to have LAPACK installed.
#Later we will make it possible if LAPACK is not present to automaticall install before compiling dendro5
if(OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()


if(MPI_COMPILE_FLAGS)
    set(COMPILE_FLAGS "${COMPILE_FLAGS} ${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
    set(LINK_FLAGS "${LINK_FLAGS} ${MPI_LINK_FLAGS}")
endif()


## options for dendro

option(USE_64BIT_INDICES "Use 64-Bit indices. Reverts to 32-bit if turned off" ON)
option(ALLTOALLV_FIX "Use K-way all to all v" OFF)
option(SPLITTER_SELECTION_FIX "Turn on Splitter Selection fix" ON)
option(DIM_2 "use the two dimentional sorting" OFF)
option(WITH_BLAS_LAPACK "build using BLAS and LAPACk" ON)
option(MANUAL_BLAS_LAPACK "configure BLAS and LAPACK Manually" OFF)
option(DENDRO_VTK_BINARY "write vtk/vtu files in binary mode " ON)
option(DENDRO_VTK_ZLIB_COMPRES "write vtk/vtu files in binary mode with zlib compression (only compatible with binary mode) " OFF)
option(BUILD_WITH_PETSC " build dendro with PETSC " ON)
option(HILBERT_ORDERING "use the Hilbert space-filling curve to order orthants" OFF)
option(BUILD_EXAMPLES "build example programs" ON)
option(PROFILING "Profile code" OFF)
set(KWAY 128 CACHE INT 128)
set(NUM_NPES_THRESHOLD 2 CACHE INT 2)

#set the build type to release by default.
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING
            "Choose the type of build, options are: Debug Release " FORCE)
endif()

if (WITH_BLAS_LAPACK)
    add_definitions(-DWITH_BLAS_LAPACK)
    if(DEFINED ENV{MKLROOT})
        find_package(LAPACK COMPONENTS MKL REQUIRED)
        set(LAPACK_LIBRARIES ${MKL_LIBRARIES})
        if (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
            set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -mkl")
        else()
            set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm")
        endif()
        message(STATUS ${LAPACK_LIBRARIES})
    elseif(MANUAL_BLAS_LAPACK)
        if ("$ENV{BLAS}" STREQUAL "")
            message("Environment Variable BLAS is not set. Please set it to BLAS directory")
        endif ()
        if ("$ENV{LAPACK}" STREQUAL "")
            message("Enviroment Variable LAPACK is note set. Please set it to LAPACK directory. ")
        endif ()
        set(LAPACKE_DIR $ENV{LAPACK}/LAPACKE)
        set(BLAS_LIBS $ENV{BLAS}/lib)
        set(LAPACK_LIBS $ENV{LAPACK}/lib)
        set(LAPACK_LINKER_FLAGS -llapacke -llapack -lblas -lgfortran -lquadmath)
        set(LAPACK_LIBRARIES ${LAPACK_LIBS}/liblapacke.a ${LAPACK_LIBS}/liblapack.a ${BLAS_LIBS}/libblas.a -static libgfortran.a libquadmath.a)
        set(LINK_FLAGS "${LINK_FLAGS} ${LAPACK_LINKER_FLAGS}")
    else ()
        find_package(BLAS REQUIRED)
        find_package(LAPACK REQUIRED)
        set(LAPACK_LINKER_FLAGS -llapacke -llapack -lblas -lgfortran -lquadmath)
        set(LAPACKE_DIR $ENV{LAPACK}/LAPACKE)
        set(LINK_FLAGS "${LINK_FLAGS} ${LAPACK_LINKER_FLAGS}")
        find_library(LAPACKE_LIB
                NAMES lapacke lapackelib liblapacke
                HINTS "/usr/lib/"
                )
        set(LAPACK_LIBRARIES ${LAPACK_LIBRARIES} ${LAPACKE_LIB})
        message(STATUS ${LAPACK_LIBRARIES})
    endif ()
endif ()

if(BUILD_WITH_PETSC)
    list (APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake-modules")
    find_package(PETSc REQUIRED)
    add_definitions(-DBUILD_WITH_PETSC)
endif()



if(DIM_2)
    add_definitions(-DDIM_2)
endif()

if(USE_64BIT_INDICES)
    add_definitions(-DUSE_64BIT_INDICES)
    #message('Configuring 64BIT indices')
endif()

if(ALLTOALLV_FIX)
    add_definitions(-DALLTOALLV_FIX)
    add_definitions(-DKWAY=${KWAY})
endif()

if(SPLITTER_SELECTION_FIX)
    add_definitions(-DSPLITTER_SELECTION_FIX)
    add_definitions(-DNUM_NPES_THRESHOLD=${NUM_NPES_THRESHOLD})
endif()

if(ALLTOALL_SPARSE)
    add_definitions(-DALLTOALL_SPARSE)
endif()

if(DENDRO_VTK_BINARY)
else()
    set(DENDRO_VTK_ZLIB_COMPRES OFF)
endif()

if(DENDRO_VTK_BINARY)
    add_definitions(-DDENDRO_VTU_BINARY)
    if(DENDRO_VTK_ZLIB_COMPRES)
        add_definitions(-DDENDRO_VTU_ZLIB)
    endif()
else()
    add_definitions(-DDENDRO_VTU_ASCII)
endif()



set(DENDRO_KT_INC include/dendro.h
                  include/dtypes.h
                  include/hcurvedata.h
                  include/KDhcurvedata_decl.h
                  include/binUtils.h
                  include/colors.h
                  include/mathUtils.h
                  include/mathUtils.tcc
                  include/octUtils.h
                  include/tnUtils.h
                  include/ompUtils.h
                  include/ompUtils.tcc
                  include/point.h
                  include/profiler.h
                  include/seqUtils.h
                  include/seqUtils.tcc
                  include/parUtils.h
                  include/parUtils.tcc
                  include/templateUtils.h
                  include/tsort.h
                  include/filterFunction.h
                  include/distTree.h
                  include/oda.h
                  include/oda.tcc
                  include/nsort.h
                  include/nsort.tcc
                  include/tsearchCmpx.h
                  include/treeNode.h
                  include/treeNode.tcc
                  include/asyncExchangeContex.h
                  include/sfcTreeLoop.h
                  include/sfcTreeLoop_matvec.h
                  include/sfcTreeLoop_matvec_io.h
                  include/meshLoop.h
                  IO/vtk/include/oct2vtk.h
                  IO/vtk/include/checkPoint.h
                  array/include/arraySlice.h
                  FEM/include/matvec.h
                  FEM/include/tensor.h
                  FEM/include/tensor_impl.h
                  FEM/include/refel.h
                  FEM/include/basis.h
                  FEM/include/lapac.h
                  FEM/include/feMat.h
                  FEM/include/feMatrix.h
                  FEM/include/feVec.h
                  FEM/include/feVector.h
                  FEM/include/interpMatrices.h
                  FEM/include/intergridTransfer.h
                  FEM/include/setDiag.h
                  FEM/include/gmgMat.h
                  )

set(DENDRO_KT_SRC src/binUtils.cpp
                  src/octUtils.cpp
                  src/parUtils.cpp
                  src/point.cpp
                  src/profiler.cpp
                  src/KDhcurvedata.cpp
                  src/KDhcurvedata_DATA.cpp
                  src/tsort.cpp
                  src/nsort.cpp
                  src/distTree.cpp
                  src/tsearchCmpx.cpp
                  src/treeNode.cpp
                  IO/vtk/src/oct2vtk.cpp
                  IO/vtk/src/checkPoint.cpp
                  src/oda.cpp
                  FEM/src/tensor.cpp
                  FEM/src/tensor_special_float.cpp
                  FEM/src/tensor_special_double.cpp
                  FEM/src/tensor_special_longlong.cpp
                  FEM/src/tensor_special_int.cpp
                  FEM/src/tensor_special_unsigned.cpp
                  FEM/src/tensor_special_char.cpp
                  FEM/src/refel.cpp
                  FEM/src/basis.cpp
                  )

if(PROFILING)
    add_definitions(-DPROFILING)
endif()

add_library(dendroKT ${DENDRO_KT_INC} ${DENDRO_KT_SRC})
target_include_directories(dendroKT PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(dendroKT PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/IO/vtk/include)
target_include_directories(dendroKT PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/IO/zlib/inc)
target_include_directories(dendroKT PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/array/include)
target_include_directories(dendroKT PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/FEM/include)
target_include_directories(dendroKT PUBLIC ${MPI_INCLUDE_PATH})
target_link_libraries(dendroKT ${MPI_LIBRARIES} m)

if(WITH_BLAS_LAPACK)
    if(MANUAL_BLAS_LAPACK)
        target_include_directories(dendroKT PUBLIC ${LAPACKE_DIR}/include)
    else()
        target_include_directories(dendroKT PUBLIC ${LAPACK_INCLUDE_DIR})
        target_include_directories(dendroKT PUBLIC ${LAPACKE_DIR}/include)
    endif()

    target_link_libraries(dendroKT ${LAPACK_LIBRARIES})

endif()

if(BUILD_WITH_PETSC)
    target_include_directories(dendroKT PUBLIC ${PETSC_INCLUDES})
    target_link_libraries(dendroKT ${PETSC_LIBRARIES})
endif()

if(HILBERT_ORDERING)
  add_definitions(-DHILBERT_ORDERING)
endif()


set(SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/test/testPartition.cpp)
add_executable(tstPartition ${SRC_FILES})
target_include_directories(tstPartition PUBLIC ${MPI_INCLUDE_PATH})
target_link_libraries(tstPartition dendroKT ${MPI_LIBRARIES} m)

set(SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/test/testConstruct.cpp)
add_executable(tstConstruct ${SRC_FILES})
target_include_directories(tstConstruct PUBLIC ${MPI_INCLUDE_PATH})
target_link_libraries(tstConstruct dendroKT ${MPI_LIBRARIES} m)


set(SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/test/testTables.cpp)
add_executable(tstTable ${SRC_FILES})
target_include_directories(tstTable PUBLIC ${MPI_INCLUDE_PATH})
target_link_libraries(tstTable dendroKT ${MPI_LIBRARIES} m)

set(SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/test/testTreeSort.cpp)
add_executable(tstTreeSort ${SRC_FILES})
target_include_directories(tstTreeSort PUBLIC ${MPI_INCLUDE_PATH})
target_link_libraries(tstTreeSort dendroKT ${MPI_LIBRARIES} m)

set(SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/test/testTreeConstruction.cpp)
add_executable(tstTreeConstruction ${SRC_FILES})
target_include_directories(tstTreeConstruction PUBLIC ${MPI_INCLUDE_PATH})
target_link_libraries(tstTreeConstruction dendroKT ${MPI_LIBRARIES} m)

set(SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/test/testTreeBalancing.cpp)
add_executable(tstTreeBalancing ${SRC_FILES})
target_include_directories(tstTreeBalancing PUBLIC ${MPI_INCLUDE_PATH})
target_link_libraries(tstTreeBalancing dendroKT ${MPI_LIBRARIES} m)

set(SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/test/testTreeVtu.cpp)
add_executable(tstTreeVtu ${SRC_FILES})
target_include_directories(tstTreeVtu PUBLIC ${MPI_INCLUDE_PATH})
target_link_libraries(tstTreeVtu dendroKT ${MPI_LIBRARIES} m)

set(SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/test/testHCurveTable.cpp)
add_executable(tstHCurveTable ${SRC_FILES})
target_include_directories(tstHCurveTable PUBLIC ${MPI_INCLUDE_PATH})
target_link_libraries(tstHCurveTable dendroKT ${MPI_LIBRARIES} m)

set(SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/test/testMovingBall.cpp)
add_executable(tstMovingBall ${SRC_FILES})
target_include_directories(tstMovingBall PUBLIC ${MPI_INCLUDE_PATH})
target_link_libraries(tstMovingBall dendroKT ${MPI_LIBRARIES} m)

set(SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/test/testTSearchCmpx.cpp)
add_executable(tstTSearchCmpx ${SRC_FILES})
target_include_directories(tstTSearchCmpx PUBLIC ${MPI_INCLUDE_PATH})
target_link_libraries(tstTSearchCmpx dendroKT ${MPI_LIBRARIES} m)

set(SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/test/testCountCGNodes.cpp)
add_executable(tstCountCGNodes ${SRC_FILES})
target_include_directories(tstCountCGNodes PUBLIC ${MPI_INCLUDE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/test/)
target_link_libraries(tstCountCGNodes dendroKT ${MPI_LIBRARIES} m)

set(SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/test/testScatterMap.cpp)
add_executable(tstScatterMap ${SRC_FILES})
target_include_directories(tstScatterMap PUBLIC ${MPI_INCLUDE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/test/)
target_link_libraries(tstScatterMap dendroKT ${MPI_LIBRARIES} m)

set(SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/test/testScatterUniq.cpp)
add_executable(tstScatterUniq ${SRC_FILES})
target_include_directories(tstScatterUniq PUBLIC ${MPI_INCLUDE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/test/)
target_link_libraries(tstScatterUniq dendroKT ${MPI_LIBRARIES} m)

set(SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/test/testODA.cpp)
add_executable(tstODA ${SRC_FILES})
target_include_directories(tstODA PUBLIC ${MPI_INCLUDE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/test/)
target_link_libraries(tstODA dendroKT ${MPI_LIBRARIES} m)

set(SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/test/testMatvec.cpp)
add_executable(tstMatvec ${SRC_FILES})
target_include_directories(tstMatvec PUBLIC ${MPI_INCLUDE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/test/)
target_link_libraries(tstMatvec dendroKT ${MPI_LIBRARIES} m)

set(SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/test/testFunction2Octree.cpp)
add_executable(tstFunction2Octree ${SRC_FILES})
target_include_directories(tstFunction2Octree PUBLIC ${MPI_INCLUDE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/test/)
target_link_libraries(tstFunction2Octree dendroKT ${MPI_LIBRARIES} m)

set(SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/test/testSFCTreeLoop.cpp)
add_executable(tstSFCTreeLoop ${SRC_FILES})
target_include_directories(tstSFCTreeLoop PUBLIC ${MPI_INCLUDE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/test/)
target_link_libraries(tstSFCTreeLoop dendroKT ${MPI_LIBRARIES} m)

set(SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/test/testKeepSiblingLeafsTogether.cpp)
add_executable(tstKeepSiblingLeafsTogether ${SRC_FILES})
target_include_directories(tstKeepSiblingLeafsTogether PUBLIC ${MPI_INCLUDE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/test/)
target_link_libraries(tstKeepSiblingLeafsTogether dendroKT ${MPI_LIBRARIES} m)

set(SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/test/testCreateRegularOctree.cpp)
add_executable(tstCreateRegularOctree ${SRC_FILES})
target_include_directories(tstCreateRegularOctree PUBLIC ${MPI_INCLUDE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/test/)
target_link_libraries(tstCreateRegularOctree dendroKT ${MPI_LIBRARIES} m)

set(SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/test/testLocGlobSz.cpp)
add_executable(tstLocGlobSz ${SRC_FILES})
target_include_directories(tstLocGlobSz PUBLIC ${MPI_INCLUDE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/test/)
target_link_libraries(tstLocGlobSz dendroKT ${MPI_LIBRARIES} m)

set(SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/test/testBillionsOfElements.cpp)
add_executable(tstBillionsOfElements ${SRC_FILES})
target_include_directories(tstBillionsOfElements PUBLIC ${MPI_INCLUDE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/test/)
target_link_libraries(tstBillionsOfElements dendroKT ${MPI_LIBRARIES} m)

set(SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/test/testMeshLoop.cpp)
add_executable(tstMeshLoop ${SRC_FILES})
target_include_directories(tstMeshLoop PUBLIC ${MPI_INCLUDE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/test/)
target_link_libraries(tstMeshLoop dendroKT ${MPI_LIBRARIES} m)

set(SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/test/testGenerateGridHierarchy.cpp)
add_executable(tstGenerateGridHierarchy ${SRC_FILES})
target_include_directories(tstGenerateGridHierarchy PUBLIC ${MPI_INCLUDE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/test/)
target_link_libraries(tstGenerateGridHierarchy dendroKT ${MPI_LIBRARIES} m)

set(SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/test/testIntergridTransfer.cpp)
add_executable(tstIntergridTransfer ${SRC_FILES})
target_include_directories(tstIntergridTransfer PUBLIC ${MPI_INCLUDE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/test/)
target_link_libraries(tstIntergridTransfer dendroKT ${MPI_LIBRARIES} m)

set(SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/test/testIntergridSubda.cpp)
add_executable(tstIntergridSubda ${SRC_FILES})
target_include_directories(tstIntergridSubda PUBLIC ${MPI_INCLUDE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/test/)
target_link_libraries(tstIntergridSubda dendroKT ${MPI_LIBRARIES} m)

set(SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/test/testRestriction.cpp)
add_executable(tstRestriction ${SRC_FILES})
target_include_directories(tstRestriction PUBLIC ${MPI_INCLUDE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/test/)
target_link_libraries(tstRestriction dendroKT ${MPI_LIBRARIES} m)

set(SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/test/testRemesh.cpp)
add_executable(tstRemesh ${SRC_FILES})
target_include_directories(tstRemesh PUBLIC ${MPI_INCLUDE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/test/)
target_link_libraries(tstRemesh dendroKT ${MPI_LIBRARIES} m)

set(SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/test/testCheckpoint.cpp)
add_executable(tstCheckpoint ${SRC_FILES})
target_include_directories(tstCheckpoint PUBLIC ${MPI_INCLUDE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/test/)
target_link_libraries(tstCheckpoint dendroKT ${MPI_LIBRARIES} m)

set(SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/test/testBoundaryIndices.cpp)
add_executable(tstBoundaryIndices ${SRC_FILES})
target_include_directories(tstBoundaryIndices PUBLIC ${MPI_INCLUDE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/test/)
target_link_libraries(tstBoundaryIndices dendroKT ${MPI_LIBRARIES} m)

set(SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/test/testRefineSubda.cpp)
add_executable(tstRefineSubda ${SRC_FILES})
target_include_directories(tstRefineSubda PUBLIC ${MPI_INCLUDE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/test/)
target_link_libraries(tstRefineSubda dendroKT ${MPI_LIBRARIES} m)

set(SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/test/testTreeNodeOwner.cpp)
add_executable(tstTreeNodeOwner ${SRC_FILES})
target_include_directories(tstTreeNodeOwner PUBLIC ${MPI_INCLUDE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/test/)
target_link_libraries(tstTreeNodeOwner dendroKT ${MPI_LIBRARIES} m)

set(SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/test/testBoundaryRefined.cpp)
add_executable(tstBoundaryRefined ${SRC_FILES})
target_include_directories(tstBoundaryRefined PUBLIC ${MPI_INCLUDE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/test/)
target_link_libraries(tstBoundaryRefined dendroKT ${MPI_LIBRARIES} m)

set(SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/test/testConsistentElements.cpp)
add_executable(tstConsistentElements ${SRC_FILES})
target_include_directories(tstConsistentElements PUBLIC ${MPI_INCLUDE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/test/)
target_link_libraries(tstConsistentElements dendroKT ${MPI_LIBRARIES} m)

set(SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/test/testProduceHangingSubda.cpp)
add_executable(tstProduceHangingSubda ${SRC_FILES})
target_include_directories(tstProduceHangingSubda PUBLIC ${MPI_INCLUDE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/test/)
target_link_libraries(tstProduceHangingSubda dendroKT ${MPI_LIBRARIES} m)

set(SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/test/testHangingBoundary.cpp)
add_executable(tstHangingBoundary ${SRC_FILES})
target_include_directories(tstHangingBoundary PUBLIC ${MPI_INCLUDE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/test/)
target_link_libraries(tstHangingBoundary dendroKT ${MPI_LIBRARIES} m)

set(SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/test/testIncidentGhostPoints.cpp
              ${CMAKE_CURRENT_SOURCE_DIR}/bench/src/genChannelPoints.cpp)
add_executable(tstIncidentGhostPoints ${SRC_FILES})
target_include_directories(tstIncidentGhostPoints PUBLIC ${MPI_INCLUDE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/test/ ${CMAKE_CURRENT_SOURCE_DIR}/bench/include)
target_link_libraries(tstIncidentGhostPoints dendroKT ${MPI_LIBRARIES} m)

set(SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/FEM/examples/include/matrix.h
              ${CMAKE_CURRENT_SOURCE_DIR}/FEM/examples/src/matrix.cpp)
add_executable(tstMatAssembly ${SRC_FILES})
set_property(TARGET tstMatAssembly PROPERTY RUNTIME_OUTPUT_DIRECTORY Examples/FEM)
target_include_directories(tstMatAssembly PUBLIC ${MPI_INCLUDE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/FEM/examples/include)
target_link_libraries(tstMatAssembly dendroKT ${MPI_LIBRARIES} m)


## tsort_bench (./tsortBench)
## -----------
set(SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/bench/include/tsort_bench.h ${CMAKE_CURRENT_SOURCE_DIR}/bench/src/tsort_bench.cpp)
add_executable(tsortBench ${SRC_FILES})
target_include_directories(tsortBench PUBLIC ${MPI_INCLUDE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/bench/include)
target_link_libraries(tsortBench dendroKT ${MPI_LIBRARIES} m)

## matvec_bench (./matvecBench)
## ------------
set(SRC_FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/bench/include/matvec_bench.h
  ${CMAKE_CURRENT_SOURCE_DIR}/bench/src/matvec_bench.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/FEM/examples/include/heatMat.h
  ${CMAKE_CURRENT_SOURCE_DIR}/FEM/examples/include/heatVec.h
  ${CMAKE_CURRENT_SOURCE_DIR}/FEM/examples/src/heatVec.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/FEM/examples/src/heatMat.cpp)
add_executable(matvecBench ${SRC_FILES})
target_include_directories(matvecBench PUBLIC ${MPI_INCLUDE_PATH}
                                              ${CMAKE_CURRENT_SOURCE_DIR}/bench/include
                                              ${CMAKE_CURRENT_SOURCE_DIR}/FEM/examples/include)
target_link_libraries(matvecBench dendroKT ${MPI_LIBRARIES} m)

## matvec_bench_adaptive (./matvecBenchAdaptive)
## ---------------------
set(SRC_FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/bench/include/matvec_bench.h
  ${CMAKE_CURRENT_SOURCE_DIR}/bench/src/matvec_bench_adaptive.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/FEM/examples/include/heatMat.h
  ${CMAKE_CURRENT_SOURCE_DIR}/FEM/examples/include/heatVec.h
  ${CMAKE_CURRENT_SOURCE_DIR}/FEM/examples/src/heatVec.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/FEM/examples/src/heatMat.cpp)
add_executable(matvecBenchAdaptive ${SRC_FILES})
target_include_directories(matvecBenchAdaptive PUBLIC ${MPI_INCLUDE_PATH}
                                              ${CMAKE_CURRENT_SOURCE_DIR}/bench/include
                                              ${CMAKE_CURRENT_SOURCE_DIR}/FEM/examples/include)
target_link_libraries(matvecBenchAdaptive dendroKT ${MPI_LIBRARIES} m)


## matvec_bench_channel (./matvecBenchChannel)
## ---------------------
set(SRC_FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/bench/include/matvec_bench.h
  ${CMAKE_CURRENT_SOURCE_DIR}/bench/src/matvec_bench_channel.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bench/include/genChannelPoints.h
  ${CMAKE_CURRENT_SOURCE_DIR}/bench/src/genChannelPoints.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/FEM/examples/include/poissonMat.h
  ${CMAKE_CURRENT_SOURCE_DIR}/FEM/examples/include/poissonVec.h
  ${CMAKE_CURRENT_SOURCE_DIR}/FEM/examples/src/poissonVec.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/FEM/examples/src/poissonMat.cpp)
add_executable(matvecBenchChannel ${SRC_FILES})
target_include_directories(matvecBenchChannel PUBLIC ${MPI_INCLUDE_PATH}
                                              ${CMAKE_CURRENT_SOURCE_DIR}/bench/include
                                              ${CMAKE_CURRENT_SOURCE_DIR}/FEM/examples/include)
target_link_libraries(matvecBenchChannel dendroKT ${MPI_LIBRARIES} m)



## Examples
if(BUILD_EXAMPLES)
  ## Heat equation.
  set(SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/FEM/examples/include/heatMat.h
                ${CMAKE_CURRENT_SOURCE_DIR}/FEM/examples/include/heatVec.h
                ${CMAKE_CURRENT_SOURCE_DIR}/FEM/examples/src/heatVec.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/FEM/examples/src/heatMat.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/FEM/examples/src/heatEq.cpp)
  add_executable(exHeatEq ${SRC_FILES})
  set_property(TARGET exHeatEq PROPERTY RUNTIME_OUTPUT_DIRECTORY Examples/FEM)
  target_include_directories(exHeatEq PUBLIC ${MPI_INCLUDE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/FEM/examples/include)
  target_link_libraries(exHeatEq dendroKT ${MPI_LIBRARIES} m)

  ## Poisson equation.
  set(SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/FEM/examples/include/poissonMat.h
                ${CMAKE_CURRENT_SOURCE_DIR}/FEM/examples/include/poissonVec.h
                ${CMAKE_CURRENT_SOURCE_DIR}/FEM/examples/src/poissonVec.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/FEM/examples/src/poissonMat.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/FEM/examples/src/poissonEq.cpp)
  add_executable(exPoissonEq ${SRC_FILES})
  set_property(TARGET exPoissonEq PROPERTY RUNTIME_OUTPUT_DIRECTORY Examples/FEM)
  target_include_directories(exPoissonEq PUBLIC ${MPI_INCLUDE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/FEM/examples/include)
  target_link_libraries(exPoissonEq dendroKT ${MPI_LIBRARIES} m)

  ## Poisson equation with gmg.
  set(SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/FEM/examples/include/poissonMat.h
                ${CMAKE_CURRENT_SOURCE_DIR}/FEM/examples/include/poissonVec.h
                ${CMAKE_CURRENT_SOURCE_DIR}/FEM/examples/src/poissonVec.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/FEM/examples/src/poissonMat.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/FEM/examples/src/poissonGMG.cpp)
  add_executable(exPoissonGMG ${SRC_FILES})
  set_property(TARGET exPoissonGMG PROPERTY RUNTIME_OUTPUT_DIRECTORY Examples/FEM)
  target_include_directories(exPoissonGMG PUBLIC ${MPI_INCLUDE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/FEM/examples/include)
  target_link_libraries(exPoissonGMG dendroKT ${MPI_LIBRARIES} m)

  ## Test the poisson elemental operators are consistent.
  set(SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/FEM/examples/include/poissonMat.h
                ${CMAKE_CURRENT_SOURCE_DIR}/FEM/examples/src/poissonMat.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/test/testPoissonOperator.cpp)
  add_executable(tstPoissonOperator ${SRC_FILES})
  set_property(TARGET tstPoissonOperator PROPERTY RUNTIME_OUTPUT_DIRECTORY Examples/FEM)
  target_include_directories(tstPoissonOperator PUBLIC ${MPI_INCLUDE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/FEM/examples/include)
  target_link_libraries(tstPoissonOperator dendroKT ${MPI_LIBRARIES} m)

endif(BUILD_EXAMPLES)


