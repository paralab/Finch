set(EXAMPLE_SRC
        src/main.cpp
        )

add_executable(adap  ${EXAMPLE_SRC})
target_include_directories(adap PUBLIC include)
target_link_libraries(adap dendriteKT dendroKT ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
