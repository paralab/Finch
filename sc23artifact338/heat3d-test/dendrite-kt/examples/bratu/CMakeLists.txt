set(EXAMPLE_INC
        include/BTEquation.h
        include/BTNodeData.h
        include/BTInputData.h
        )

set(EXAMPLE_SRC
        src/bratu.cpp
        )

add_executable(bratu  ${EXAMPLE_SRC} ${EXAMPLE_INC})
target_include_directories(bratu PUBLIC include)
target_link_libraries(bratu dendriteKT dendroKT  ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
