set(EXAMPLE_INC
        include/SSHTEquation.h
        include/SSHTNodeData.h
        include/SSHTInputData.h
        )

set(EXAMPLE_SRC
        src/SSHT.cpp
        )

add_executable(steady-heat  ${EXAMPLE_SRC} ${EXAMPLE_INC})
target_include_directories(steady-heat PUBLIC include)
target_link_libraries(steady-heat dendriteKT dendroKT  ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
