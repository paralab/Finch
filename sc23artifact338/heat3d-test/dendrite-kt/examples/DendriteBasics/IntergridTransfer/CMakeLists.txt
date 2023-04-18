
set(EXAMPLE_INC
    include/DARefine.h
    include/DARefineValues.h
        )


set(EXAMPLE_SRC
        src/main
        )

add_executable(intergrid ${EXAMPLE_SRC} ${EXAMPLE_INC})
target_include_directories(intergrid PUBLIC include)
target_link_libraries(intergrid dendriteKT dendroKT  ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
