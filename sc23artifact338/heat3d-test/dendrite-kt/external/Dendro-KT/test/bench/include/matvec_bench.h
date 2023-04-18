/**
 * @brief: Benchmark program to time 4D matvec:
 * - [Adaptive grid] treesort octree construction, balancing, create oda (cg node identification + scatter map),
 * - [Regular grid] matvec: communication, top-down pass, bottom-up pass, leaf-node elemental matvec.
 * @date: 04/09/2019
 * @author: Milinda Fernando, School of Computing University of Utah. 
 * @author: Masado Ishii, School of Computing University of Utah. 
 * 
 * @note: Based on FEM/examples/src/heatEq.cpp and bench/src/tsort_bench.cpp
*/


#ifndef DENDRO_KT_MATVEC_BENCH_H
#define DENDRO_KT_MATVEC_BENCH_H

#include <vector>
#include <assert.h>
#include <mpi.h>
#include "profiler.h"
#include <iostream>


    namespace bench
    {
        extern profiler_t t_adaptive_tsort;
        extern profiler_t t_adaptive_tconstr;
        extern profiler_t t_adaptive_tbal;
        extern profiler_t t_adaptive_oda;

        // For now, I will inject timing into the implementation.
        extern profiler_t t_ghostexchange;
        extern profiler_t t_topdown;
        extern profiler_t t_bottomup;
        extern profiler_t t_treeinterior;
        extern profiler_t t_elemental;
        extern profiler_t t_matvec;
        
        /**@breif reset all the extern counters defined here.. */
        void resetAllTimers();


        /**
         * @brief: performs the bench kernel node sort, tree construction and balancing. 
         * @param [in] numPts: number of points/elements per core, 
         * @param [in] numIter: number of iterations
         * @param [in] comm: MPI communicator
         * 
        */
        void bench_kernel(unsigned int numPts, unsigned int numIter = 10, MPI_Comm comm=MPI_COMM_WORLD);

        void dump_profile_info(std::ostream& fout, profiler_t* timers, char ** names, unsigned int n ,MPI_Comm comm=MPI_COMM_WORLD);

    } // end of namespace bench. 


#endif



