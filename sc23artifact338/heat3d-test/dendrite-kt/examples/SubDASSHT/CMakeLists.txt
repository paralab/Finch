set(EXAMPLE_INC
        include/SSHTEquation.h
        include/SSHTNodeData.h
        include/SSHTAnalytic.h
        )

set(EXAMPLE_SRC
        src/SSHT.cpp
        )

add_executable(subDA-steady-heat  ${EXAMPLE_SRC} ${EXAMPLE_INC})
target_include_directories(subDA-steady-heat PUBLIC include)
target_link_libraries(subDA-steady-heat dendriteKT dendroKT  ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
