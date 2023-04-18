set(EXAMPLE_INC
        include/Refine.h
        )

set(EXAMPLE_SRC
        src/main.cpp
        )

add_executable(ibm-compar  ${EXAMPLE_SRC} ${EXAMPLE_INC})
target_include_directories(ibm-compar  PUBLIC include)
target_link_libraries(ibm-compar  dendriteKT dendroKT  ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
