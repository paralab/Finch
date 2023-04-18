set(EXAMPLE_SRC
        src/main.cpp
        )

add_executable(sda  ${EXAMPLE_SRC})
target_include_directories(sda PUBLIC include)
target_link_libraries(sda dendriteKT dendroKT ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
