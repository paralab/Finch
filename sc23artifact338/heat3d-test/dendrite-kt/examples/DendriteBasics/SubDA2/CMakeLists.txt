
set(EXAMPLE_INC
#        include/Loop.h
#        include/SSHTEquation.h
#        include/SSHTNodeData.h
        include/SDARefine.h
        )


set(EXAMPLE_SRC
        src/main.cpp
        )

add_executable(sda2 ${EXAMPLE_SRC} ${EXAMPLE_INC})
target_include_directories(sda2 PUBLIC include)
target_link_libraries(sda2 dendriteKT dendroKT  ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
