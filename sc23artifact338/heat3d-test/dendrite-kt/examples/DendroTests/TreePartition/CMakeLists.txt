set(EXAMPLE_SRC
        src/main.cpp
        )

add_executable(treePart  ${EXAMPLE_SRC})
target_include_directories(treePart PUBLIC include)
target_link_libraries(treePart dendriteKT dendroKT ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
