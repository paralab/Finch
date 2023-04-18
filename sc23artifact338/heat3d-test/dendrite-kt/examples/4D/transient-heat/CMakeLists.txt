set(EXAMPLE_INC
        include/TSHTEquation.h
        include/TSHTNodeData.h
        include/TSHTAnalytic.h
        )

set(EXAMPLE_SRC
        src/TSHT.cpp
        )

add_executable(tsht  ${EXAMPLE_SRC} ${EXAMPLE_INC})
target_include_directories(tsht PUBLIC include)
target_link_libraries(tsht dendriteKT dendroKT  ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
