set(EXAMPLE_INC
        include/SSHTEquation.h
        include/SSHTNodeData.h
        )

set(EXAMPLE_SRC
        src/SSHT.cpp
        )

add_executable(ssht-neumann  ${EXAMPLE_SRC} ${EXAMPLE_INC})
target_include_directories(ssht-neumann PUBLIC include)
target_link_libraries(ssht-neumann dendriteKT dendroKT  ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
