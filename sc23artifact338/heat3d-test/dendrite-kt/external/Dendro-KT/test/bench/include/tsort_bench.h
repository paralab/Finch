/**
 * @brief: Simple benchmark program to time 
 * treesort octree construction and balancing with cg node identification and 
 * building the scatter map. 
 * @date: 03/28/1019
 * @author: Milinda Fernando, School of Computing University of Utah. 
 * 
*/

#ifndef BENDRO_KT_TSORT_BENCH_H
#define BENDRO_KT_TSORT_BENCH_H

#include "treeNode.h"
#include "tsort.h"
#include "nsort.h"
#include "octUtils.h"
#include "hcurvedata.h"
#include "octUtils.h"
#include <vector>
#include <assert.h>
#include <mpi.h>
#include "profiler.h"
#include <iostream>


    namespace bench
    {
        /**@brief: begin of tree sort*/
        extern profiler_t t_sort;
        
        /**@brief:begin of 2:1 balance*/
        extern profiler_t t_con;

        /**@brief:begin of 2:1 balance*/
        extern profiler_t t_bal;
        
        /**@brief: begin to cg node compute*/
        extern profiler_t t_cg;
        
        /**@brief: scatter map begin*/
        extern profiler_t t_sm;
        
        /**@breif reset all the extern counters defined here.. */
        void resetAllTimers();


        /**
         * @brief: performs the bench kernel node sort, tree construction and balancing. 
         * @param [in] number of points per core, 
         * @param [in] numIter: number of iterations
         * @param [in] comm: MPI communicator
         * 
        */
        void bench_kernel(unsigned int numPts, unsigned int numIter,unsigned int pOrder, MPI_Comm comm=MPI_COMM_WORLD);

        void dump_profile_info(std::ostream& fout, profiler_t* timers, char ** names, unsigned int n ,MPI_Comm comm=MPI_COMM_WORLD);

    } // end of namespace bench. 


#endif



