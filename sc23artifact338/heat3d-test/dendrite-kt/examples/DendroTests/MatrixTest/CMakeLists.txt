set(EXAMPLE_INC
        include/matrix.h
        include/vector.h
        )

set(EXAMPLE_SRC
        src/main.cpp
        )

add_executable(matTest  ${EXAMPLE_INC} ${EXAMPLE_SRC})
target_include_directories(matTest PUBLIC include)
target_link_libraries(matTest dendriteKT dendroKT ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
