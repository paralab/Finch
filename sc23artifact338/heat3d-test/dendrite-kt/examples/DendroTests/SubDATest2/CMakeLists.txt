set(EXAMPLE_SRC
        src/main.cpp
        )

add_executable(sda-test-dendro  ${EXAMPLE_SRC})
target_include_directories(sda-test-dendro PUBLIC include)
target_link_libraries(sda-test-dendro  dendroKT dendriteKT ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
