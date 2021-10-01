cmake_minimum_required(VERSION 2.8)
project(poisson3d)

set(poisson3d_INC include/linear_skel.h
            include/bilinear_skel.h
            include/Genfunction.h
        )

set(poisson3d_SRC src/linear_skel.cpp
            src/bilinear_skel.cpp
            src/Genfunction.cpp
        )

set(SOURCE_FILES src/poisson3d.cpp ${poisson3d_INC} ${poisson3d_SRC})
add_executable(poisson3d ${SOURCE_FILES})
target_include_directories(poisson3d PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(poisson3d PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_include_directories(poisson3d PRIVATE ${CMAKE_SOURCE_DIR}/include/test)
target_include_directories(poisson3d PRIVATE ${CMAKE_SOURCE_DIR}/examples/include)
target_include_directories(poisson3d PRIVATE ${CMAKE_SOURCE_DIR}/FEM/include)
target_include_directories(poisson3d PRIVATE ${CMAKE_SOURCE_DIR}/ODE/include)
target_include_directories(poisson3d PRIVATE ${CMAKE_SOURCE_DIR}/LinAlg/include)
target_include_directories(poisson3d PRIVATE ${CMAKE_SOURCE_DIR}/IO/vtk/include)
target_include_directories(poisson3d PRIVATE ${CMAKE_SOURCE_DIR}/IO/zlib/inc)
target_include_directories(poisson3d PRIVATE ${MPI_INCLUDE_PATH})
target_include_directories(poisson3d PRIVATE ${GSL_INCLUDE_DIRS})
if(WITH_CUDA)
    target_include_directories(poisson3d PRIVATE ${CUDA_INCLUDE_DIRS})
endif()
target_link_libraries(poisson3d dendro5 ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)

