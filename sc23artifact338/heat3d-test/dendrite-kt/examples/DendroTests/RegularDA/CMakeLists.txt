set(EXAMPLE_SRC
        src/main.cpp
        )

add_executable(rda  ${EXAMPLE_SRC})
target_include_directories(rda PUBLIC include)
target_link_libraries(rda dendriteKT dendroKT ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
