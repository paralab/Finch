
set(EXAMPLE_INC
#        include/SSHTEquation.h
#        include/SSHTNodeData.h
        include/SDARefine.h
        )


set(EXAMPLE_SRC
        src/main.cpp
        )

add_executable(sda-geom ${EXAMPLE_SRC} ${EXAMPLE_INC})
target_include_directories(sda-geom PUBLIC include)
target_link_libraries(sda-geom dendriteKT dendroKT  ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
