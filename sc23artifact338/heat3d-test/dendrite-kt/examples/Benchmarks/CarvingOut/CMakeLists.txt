set(EXAMPLE_INC
        include/SDARefine.h
        )

set(EXAMPLE_SRC
        src/main.cpp
        )

add_executable(bench-carved  ${EXAMPLE_SRC} ${EXAMPLE_INC})
target_include_directories(bench-carved  PUBLIC include)
target_link_libraries(bench-carved  dendriteKT dendroKT  ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
