set(EXAMPLE_INC
        include/SSHTNodeData.h
        include/SSHTNodeData.h
        include/CalcError.h
        include/IMGALoop.h
        )

set(EXAMPLE_SRC
        src/main.cpp
        )

add_executable(imga-test  ${EXAMPLE_SRC} ${EXAMPLE_INC})
target_include_directories(imga-test PUBLIC include)
target_link_libraries(imga-test dendriteKT dendroKT  ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
