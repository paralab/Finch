
set(EXAMPLE_INC
        include/SDARefine.h
        )


set(EXAMPLE_SRC
        src/checkpoint.cpp
        )

add_executable(checkpoint  ${EXAMPLE_SRC} ${EXAMPLE_INC} )
target_include_directories(checkpoint PUBLIC include)
target_link_libraries(checkpoint dendriteKT dendroKT  ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
