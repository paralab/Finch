set(EXAMPLE_INC
        include/MomentumEquation.h
        include/NSNodeData.h
        include/NSRefine.h
        include/NSUtils.h
        include/NSInputData.h
        include/NSBoundaryConditions.h
        include/PressurePoissonEquation.h
        include/VelocityUpdateEquation.h
        include/NSParams.h
        )

set(EXAMPLE_SRC
        src/main.cpp
        )

add_executable(ns-proj  ${EXAMPLE_SRC} ${EXAMPLE_INC})
target_include_directories(ns-proj PUBLIC include)
target_link_libraries(ns-proj dendriteKT dendroKT  ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
