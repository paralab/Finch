set(EXAMPLE_INC
#        include/matrix.h
#        include/vector.h
        )

set(EXAMPLE_SRC
        src/main.cpp
        )

add_executable(intergrid-dendro  ${EXAMPLE_INC} ${EXAMPLE_SRC})
target_include_directories(intergrid-dendro PUBLIC include)
target_link_libraries(intergrid-dendro dendriteKT dendroKT ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
