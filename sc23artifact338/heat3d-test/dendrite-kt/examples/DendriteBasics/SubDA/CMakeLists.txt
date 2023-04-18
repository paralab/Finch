
set(EXAMPLE_INC
        include/SSHTEquation.h
        include/SSHTNodeData.h
        include/SDARefine.h)


set(EXAMPLE_SRC
        src/main.cpp
        )

add_executable(sda-test ${EXAMPLE_SRC} ${EXAMPLE_INC})
target_include_directories(sda-test PUBLIC include)
target_link_libraries(sda-test dendriteKT dendroKT  ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
