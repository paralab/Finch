set(EXAMPLE_INC
        include/SDARefine.h
        )

set(EXAMPLE_SRC
        src/main.cpp
        )

add_executable(bench-channel  ${EXAMPLE_SRC} ${EXAMPLE_INC})
target_include_directories(bench-channel  PUBLIC include)
target_link_libraries(bench-channel dendriteKT dendroKT  ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
