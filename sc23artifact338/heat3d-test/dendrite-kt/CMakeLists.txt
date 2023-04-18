## CMAKE for Dendro-KT
cmake_minimum_required(VERSION 2.8)
project(dendriteKT)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -ftemplate-depth-2000")
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake-modules")  # for FindPETSc
find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)


# For now we just make it compulsory to have LAPACK installed.
#Later we will make it possible if LAPACK is not present to automaticall install before compiling dendro5
if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif ()


if (MPI_COMPILE_FLAGS)
    set(COMPILE_FLAGS "${COMPILE_FLAGS} ${MPI_COMPILE_FLAGS}")
endif ()

if (MPI_LINK_FLAGS)
    set(LINK_FLAGS "${LINK_FLAGS} ${MPI_LINK_FLAGS}")
endif ()


## options for dendro

option(USE_64BIT_INDICES "Use 64-Bit indices. Reverts to 32-bit if turned off" ON)
option(ALLTOALLV_FIX "Use K-way all to all v" OFF)
option(SPLITTER_SELECTION_FIX "Turn on Splitter Selection fix" ON)
option(DIM_2 "use the two dimentional sorting" OFF)
option(WITH_BLAS_LAPACK "build using BLAS and LAPACk" ON)
option(MANUAL_BLAS_LAPACK "configure BLAS and LAPACK Manually" OFF)
option(DENDRO_VTK_BINARY "write vtk/vtu files in binary mode " ON)
option(DENDRITE_VTU_ASCII "write vtk/vtu files in ASCII mode " OFF)
option(DENDRO_VTK_ZLIB_COMPRES "write vtk/vtu files in binary mode with zlib compression (only compatible with binary mode) " OFF)
option(BUILD_WITH_PETSC " build dendro with PETSC " ON)
option(HILBERT_ORDERING "use the Hilbert space-filling curve to order orthants" OFF)
option(BUILD_EXAMPLES "build example programs" ON)

## Dimension
option(ENABLE_4D "enable 4D computation" OFF)
option(ENABLE_2D "enable 2D computation" OFF)
option(ENABLE_3D "enable 3D computation" ON)

# Tensor computation
option(TENSOR "Use Tensor Operation" OFF)
option(PROFILING "Enable profiling" OFF)

# imga
option(IBM "Enable IBM Functions" OFF)

set(KWAY 128 CACHE INT 128)
set(NUM_NPES_THRESHOLD 2 CACHE INT 2)

#set the build type to release by default.
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Debug" CACHE STRING
            "Choose the type of build, options are: Debug Release " FORCE)
endif ()

if (WITH_BLAS_LAPACK)
    add_definitions(-DWITH_BLAS_LAPACK)
    if (DEFINED ENV{MKLROOT})
        message("MKL Activated")
        find_package(LAPACK COMPONENTS MKL REQUIRED)
        set(LAPACK_LIBRARIES ${MKL_LIBRARIES})
        if (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
            set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -mkl")
        else ()
            set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm")
        endif ()
        message(STATUS ${LAPACK_LIBRARIES})
    elseif (MANUAL_BLAS_LAPACK)
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

if (BUILD_WITH_PETSC)
    find_package(PETSc REQUIRED)
    add_definitions(-DBUILD_WITH_PETSC)
endif ()

if (TENSOR)
    add_definitions(-DTENSOR)
    message("Enabling Tensor operation")
endif ()
if (DIM_2)
    add_definitions(-DDIM_2)
endif ()

if (PROFILING)
    message("Enabled kT profiling")
    add_definitions(-DPROFILING)
endif ()
if (USE_64BIT_INDICES)
    add_definitions(-DUSE_64BIT_INDICES)
    #message('Configuring 64BIT indices')
endif ()

if (ALLTOALLV_FIX)
    add_definitions(-DALLTOALLV_FIX)
    add_definitions(-DKWAY=${KWAY})
endif ()

if (SPLITTER_SELECTION_FIX)
    add_definitions(-DSPLITTER_SELECTION_FIX)
    add_definitions(-DNUM_NPES_THRESHOLD=${NUM_NPES_THRESHOLD})
endif ()

if (ALLTOALL_SPARSE)
    add_definitions(-DALLTOALL_SPARSE)
endif ()

if (DENDRO_VTK_BINARY)
else ()
    set(DENDRO_VTK_ZLIB_COMPRES OFF)
endif ()

if (DENDRITE_VTU_ASCII)
    add_definitions(-DDENDRITE_VTU_ASCII)
    message("Writing in ASCII format")
else ()
    message("Writing in Binary format")
endif ()


if (DENDRO_VTK_BINARY)
    add_definitions(-DDENDRO_VTU_BINARY)
    if (DENDRO_VTK_ZLIB_COMPRES)
        add_definitions(-DDENDRO_VTU_ZLIB)
    endif ()
else ()
    add_definitions(-DDENDRO_VTU_ASCII)
endif ()

if (ENABLE_4D)
    add_definitions(-DENABLE_4D)
    message("Enabling 4D computation")
    set(ENABLE_4D,ON)
    set(ENABLE_3D,OFF)
elseif (ENABLE_2D)
    add_definitions(-DENABLE_2D)
    message("Enabling 2D computation")
    set(ENABLE_2D,ON)
    set(ENABLE_3D,OFF)
elseif (ENABLE_3D)
    add_definitions(-DENABLE_3D)
    set(ENABLE_3D,ON)
    message("Enabling 3D computation")
endif ()

if (IBM)
    add_definitions(-DIBM)
    message("Enabling IMGA Computation")
endif ()

set(EXTERNAL_TALY_LOCATION ${CMAKE_CURRENT_SOURCE_DIR}/external/talylite)
add_subdirectory(${EXTERNAL_TALY_LOCATION})
set(TALYFEM_BUILD_TUTORIALS OFF)
set(TALYFEM_BUILD_TESTS OFF)

if (HILBERT_ORDERING)
    add_definitions(-DHILBERT_ORDERING)
endif ()


#external subdirectory
add_subdirectory(external/Dendro-KT)
add_subdirectory(external/ray_tracing)
set(DendroTests OFF)


set(DENDRITE_INC
        include/DendriteUtils.h
        include/DataTypes.h

        include/TalyDendroSync.h
        include/TalyMesh.h
        include/TalyVec.h
        include/TalyMat.h
        include/TalyEquation.h
        include/MatVecCommon.h
        include/OctToPhysical.h
        include/NodeAndValues.h
        include/NodeAndValues.h

        #VMS
        include/VMSparams.h
        #Time
        include/TimeInfo.h

        #Solvers
        include/PETSc/Solver/Solver.h
        include/PETSc/BoundaryConditions.h
        include/PETSc/Solver/LinearSolver.h
        include/PETSc/Solver/NonLinearSolver.h
        include/PETSc/PetscUtils.h
        include/PETSc/VecInfo.h
        include/PETSc/VecBounds.h
        include/PETSc/VecUtils.h

        #IO files
        include/IO/VTU.h
        include/PETSc/IO/petscVTU.h
        include/PointCloud/PointData.h

        #Traversal
        include/Traversal/Traversal.h
        include/Traversal/Refinement.h
        include/Traversal/Analytic.h
        include/Traversal/DomainBounds.h
        include/Traversal/SurfaceLoop.h

        #Boundary
        include/Boundary/DomainBoundary.h
        include/Boundary/BoundaryDataTypes.h

        #Checkpointing
        include/Checkpoint/Checkpointer.h
        )
set(DENDRITE_SRC
        src/DendriteUtils.cpp
        src/OctToPhysical.cpp

        src/PETSc/Solver/Solver.cpp
        src/PETSc/Solver/LinearSolver.cpp
        src/PETSc/Solver/NonLinearSolver.cpp
        src/PETSc/BoundaryConditions.cpp
        src/PETSc/VecInfo.cpp

        #IO files
        src/IO/VTU.cpp
        src/PETSc/IO/petscVTU.cpp

        src/Traversal/Traversal.cpp
        src/Traversal/Refinement.cpp
        src/Traversal/Analytic.cpp
        src/Traversal/DomainBounds.cpp
        src/Traversal/SurfaceLoop.cpp

        #Boundary
        src/Boundary/DomainBoundary.cpp

        #Checkpointing
        src/Checkpoint/Checkpointer.cpp
        )

set(SUBDA_INC
        include/SubDA/Voxel.h
        include/SubDA/SubDomain.h
        include/Boundary/SubDomainBoundary.h)


set(SUBDA_SRC
        src/SubDA/Voxel.cpp
        src/SubDA/SubDomain.cpp
        src/Boundary/SubDomainBoundary.cpp
        )

set(GEOMETRY_INC
        include/Geometry/STL.h
        include/Geometry/MSH.h
        include/Geometry/Geometry.h)
set(GEOMETRY_SRC
        src/Geometry/STL.cpp
        src/Geometry/MSH.cpp
        src/Geometry/Geometry.cpp)

if (TENSOR)
    set(TENSOR_INC
            include/Basis/Basis.h
            include/Basis/MatVec.h
            include/Basis/Vec.h
            include/Basis/Mat.h
            )

endif ()

if (IBM)
        set(IBM_INC
                include/IMGA/IMGAUtils.h
                include/IMGA/Geometry/Geometry.h
                include/IMGA/Geometry/GeomType.h
                include/IMGA/IMGA.h
                include/IMGA/IMGADataTypes.h
                include/IMGA/Marker.h
                include/IMGA/IMGATraversal.h
                include/IMGA/SurfaceValues.h
                include/IMGA/IMGASolverUtils.h
                include/IMGA/IMGAInterpolate.h
                include/IMGA/MovingIMGA.h
                include/IMGA/MathOperation.h
                include/IMGA/SixDof.h
                include/IMGA/IO.h
                )

        set(IBM_SRC
                src/IMGA/Geometry/Geometry.cpp
                src/IMGA/IMGA.cpp
                src/IMGA/MathOperation.cpp
                )

endif ()

if(PROFILING)
    set(PROFILING_INC
            include/Profiling/GlobalProfiler.h
            )
endif()
add_library(dendriteKT ${DENDRITE_INC} ${DENDRITE_SRC} ${IBM_INC} ${IBM_SRC} ${TENSOR_INC} ${SUBDA_INC} ${SUBDA_SRC}
        ${GEOMETRY_INC} ${GEOMETRY_SRC} ${PROFILING_INC})
target_include_directories(dendriteKT PUBLIC include)

target_link_libraries(dendriteKT dendroKT talyfem raytracer ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
option(DENDRITEkT_BUILD_EXAMPLES "Build example projects" ON)
if (DENDRITEkT_BUILD_EXAMPLES)
    if (ENABLE_4D)
        add_subdirectory(examples/4D/transient-heat)
        add_subdirectory(examples/SSHT)
    else ()
        add_subdirectory(examples/SSHT)
        add_subdirectory(examples/SSHT-neumann)
        add_subdirectory(examples/DendroTests/RegularDA)
        add_subdirectory(examples/DendroTests/SubDA)
        add_subdirectory(examples/DendroTests/AdaptiveMesh)
        add_subdirectory(examples/DendriteBasics/SubDA)
        add_subdirectory(examples/DendroTests/MatrixTest)
        add_subdirectory(examples/DendroTests/IntergridTransfer)
        add_subdirectory(examples/Benchmarks/CarvingOut)
        add_subdirectory(examples/TSHT)
        add_subdirectory(examples/PredatorPrey)
        add_subdirectory(examples/bratu)
        add_subdirectory(examples/Navier-Stokes-PSPG)
        add_subdirectory(examples/Navier-Stokes-Projection)
        add_subdirectory(examples/DendriteBasics/PointCloud)
        add_subdirectory(examples/DendriteBasics/Checkpointing)
        add_subdirectory(examples/DendriteBasics/SurfaceLoop)
        add_subdirectory(examples/DendriteBasics/SubDA2)
        add_subdirectory(examples/DendriteBasics/IntergridTransfer)
        add_subdirectory(examples/DendroTests/SubDATest2)
        add_subdirectory(examples/Benchmarks/IBMComparison)
        add_subdirectory(examples/Benchmarks/Channel)
        if (IBM)
            add_subdirectory(examples/DendroTests/TreePartition)
            add_subdirectory(examples/IMGA/Basic)
        endif ()
    endif ()
endif ()

if (ENABLE_3D)
    add_subdirectory(examples/DendriteBasics/Geometry)
endif ()
