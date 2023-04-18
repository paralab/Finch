
set(EXAMPLE_INC
        include/Loop.h
#        include/SSHTEquation.h
#        include/SSHTNodeData.h
#        include/SDARefine.h
        )


set(EXAMPLE_SRC
        src/main.cpp
        )

add_executable(sloop ${EXAMPLE_SRC} ${EXAMPLE_INC})
target_include_directories(sloop PUBLIC include)
target_link_libraries(sloop dendriteKT dendroKT  ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
