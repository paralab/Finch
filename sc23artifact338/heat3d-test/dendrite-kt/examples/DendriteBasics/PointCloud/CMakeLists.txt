

set(EXAMPLE_SRC
        src/main.cpp
        )

add_executable(pcloud ${EXAMPLE_SRC})
target_include_directories(pcloud PUBLIC include)
target_link_libraries(pcloud dendriteKT dendroKT  ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
