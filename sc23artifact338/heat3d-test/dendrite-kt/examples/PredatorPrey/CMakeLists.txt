set(EXAMPLE_INC
        include/PPEquation.h
        include/PPNodeData.h
        )

set(EXAMPLE_SRC
        src/PP.cpp
        )

add_executable(predator-prey  ${EXAMPLE_SRC} ${EXAMPLE_INC})
target_include_directories(predator-prey PUBLIC include)
target_link_libraries(predator-prey dendriteKT dendroKT  ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
