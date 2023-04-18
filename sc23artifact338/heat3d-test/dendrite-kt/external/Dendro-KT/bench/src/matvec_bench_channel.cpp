/**
 * @brief: Benchmark program to time 4D matvec:
 * - treesort octree construction, balancing, create oda (cg node identification + scatter map),
 * - matvec: communication, top-down pass, bottom-up pass, leaf-node elemental matvec.
 * - In matvec_bench_adaptive, an adaptive grid is used for the matvec.
 * @date: 04/10/2019
 * @author: Milinda Fernando, School of Computing University of Utah. 
 * @author: Masado Ishii, School of Computing University of Utah. 
 * 
 * @note: Based on FEM/examples/src/poissonEq.cpp and bench/src/tsort_bench.cpp
*/

#include "matvec_bench.h"
#include "genChannelPoints.h"

#include "treeNode.h"
#include "tsort.h"
#include "nsort.h"
#include "octUtils.h"
#include "hcurvedata.h"
#include "filterFunction.h"
#include "distTree.h"

#include "matvec.h"
#include "feMatrix.h"

#include "poissonMat.h"
#include "poissonVec.h"

#include <cstring>


// ======================
// Signal handler
// https://stackoverflow.com/a/77336
// ======================

#include <stdio.h>
#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>

void handler(int sig) {
  void *array[20];
  size_t size;

  // get void*'s for all entries on the stack
  size = backtrace(array, 20);

  // print out all the frames to stderr
  fprintf(stderr, "Error: signal %d:\n", sig);
  backtrace_symbols_fd(array, size, STDERR_FILENO);
  exit(1);
}
// ======================

namespace bench
{
    profiler_t t_adaptive_tsort;
    profiler_t t_adaptive_tconstr;
    profiler_t t_adaptive_tbal;
    profiler_t t_adaptive_oda;

    profiler_t t_ghostexchange;
    profiler_t t_topdown;
    profiler_t t_bottomup;
    profiler_t t_treeinterior;
    profiler_t t_elemental;
    profiler_t t_matvec;

    struct ReportSizes
    {
      unsigned int b1_treeSortSz;
      unsigned int b1_treeConstructionSz;
      unsigned int b1_treeBalancingSz;
      unsigned int b2_treeMatvecSz;

      unsigned int b1_globNodeSz;    // Only used globally.
      unsigned int b2_globNodeSz;    // Only used globally.
    } gRptSz, gDistRptSz;


    void resetAllTimers()
    {
        t_adaptive_tsort.clear();
        t_adaptive_tconstr.clear();
        t_adaptive_tbal.clear();
        t_adaptive_oda.clear();

        t_ghostexchange.clear();
        t_topdown.clear();
        t_bottomup.clear();
        t_treeinterior.clear();
        t_elemental.clear();
        t_matvec.clear();
        
    }



    template <unsigned int dim>
    void bench_kernel(unsigned int numPts, unsigned int numTreeRuns, unsigned int numWarmup, unsigned int numRuns, unsigned int eleOrder, int lengthPower2, MPI_Comm comm)
    {
        // numWarmup affects number of matVec warmup runs.
        // numRuns affects number of matVec runs.
        // The construction/balancing example executes 1 warmup run and 1 live run.

        int rank, npes;
        MPI_Comm_rank(comm,&rank);
        MPI_Comm_size(comm,&npes);

        resetAllTimers();
        using T = unsigned int;
        using TreeNode = ot::TreeNode<T,dim>;
        using TNP = ot::TNPoint<T,dim>;
        using SFC_Tree = ot::SFC_Tree<T,dim>;

        const unsigned int maxPtsPerRegion = 1;
        const double loadFlexibility = 0.1;

        /// const int lengthPower2 = 4;  // 16 = 2^4;

        const ibm::DomainDecider boxDecider = getBoxDecider<dim>(lengthPower2);

        std::vector<ot::TreeNode<unsigned int, dim>> points = getChannelPoints<dim>(
            numPts, lengthPower2, comm);

        // Remove duplicates that could force leafs to m_uiMaxDepth.
        ot::SFC_Tree<T,dim>::distRemoveDuplicates(points, loadFlexibility, false, comm);

        ot::DistTree<unsigned int, dim> dtree;

        // 1. Benchmark construction/balancing/ODA on adaptive grid.
        {
            std::vector<TreeNode> tree;

            // Warmup run for adaptive grid.
            std::vector<TreeNode> points_copy = points;
            ot::SFC_Tree<T,dim>::distTreeBalancingWithFilter(boxDecider, points_copy, tree, maxPtsPerRegion, loadFlexibility, comm);

            //
            // Benchmark the adaptive grid example.
            //

            // Time sorting.
            for (int ii = 0; ii < numTreeRuns; ii++)
            {
              points_copy = points;
              tree.clear();
              t_adaptive_tsort.start();
              ot::SFC_Tree<T,dim>::distTreeSort(points_copy, loadFlexibility, comm);
              t_adaptive_tsort.stop();
            }
            gRptSz.b1_treeSortSz = points_copy.size();

            // Time construction.
            for (int ii = 0; ii < numTreeRuns; ii++)
            {
              points_copy = points;
              tree.clear();
              t_adaptive_tconstr.start();
              ot::SFC_Tree<T,dim>::distTreeConstructionWithFilter(boxDecider, points_copy, tree, maxPtsPerRegion, loadFlexibility, comm);
              t_adaptive_tconstr.stop();
            }
            gRptSz.b1_treeConstructionSz = tree.size();

            // Time balanced construction.
            for (int ii = 0; ii < numTreeRuns; ii++)
            {
              points_copy = points;
              tree.clear();
              t_adaptive_tbal.start();
              ot::SFC_Tree<T,dim>::distTreeBalancingWithFilter(boxDecider, points_copy, tree, maxPtsPerRegion, loadFlexibility, comm);
              t_adaptive_tbal.stop();
            }
            gRptSz.b1_treeBalancingSz = tree.size();

            dtree = ot::DistTree<unsigned int, dim>(tree, comm);
            dtree.filterTree(boxDecider);

            // Generate DA from balanced tree.
            for (int ii = 0; ii < numTreeRuns; ii++)
            {
              t_adaptive_oda.start();
              ot::DA<dim> oda(dtree, comm, eleOrder, numPts, loadFlexibility);
              t_adaptive_oda.stop();
              gDistRptSz.b1_globNodeSz = oda.getGlobalNodeSz();
            }
        }

        // 2. Benchmark matvec on regular grid.
        {
            // I've removed the right-hand side and iterative solving.
            // In this pass all we do is execute matvec in a loop.

            // Construct an adaptive grid based on a Gaussian point cloud.
            std::vector<TreeNode> tree = dtree.getTreePartFiltered();
            gRptSz.b2_treeMatvecSz = tree.size();

            // DA based on adaptive grid.
            ot::DA<dim> *octDA = new ot::DA<dim>(dtree, comm, eleOrder, numPts, loadFlexibility);
            gDistRptSz.b2_globNodeSz = octDA->getGlobalNodeSz();

            const unsigned int DOF = 1;   // matvec only supports dof==1 right now.

            std::vector<double> uSolVec, fVec, mfVec, dummyVec;
            octDA->createVector(uSolVec,false,false,DOF);
            /// octDA->createVector(fVec,false,false,DOF);
            /// octDA->createVector(mfVec,false,false,DOF);
            octDA->createVector(dummyVec,false,false,DOF);
            double *uSolVecPtr=&(*(uSolVec.begin()));
            /// double *fVecPtr=&(*(fVec.begin()));
            /// double *mfVecPtr=&(*(mfVec.begin()));
            double *dummyVecPtr=&(*(dummyVec.begin()));

            Point<dim> domain_min(-0.5,-0.5,-0.5);
            Point<dim> domain_max(0.5,0.5,0.5);

            PoissonEq::PoissonMat<dim> myPoissonMat(octDA, &tree,DOF);
            myPoissonMat.setProblemDimensions(domain_min,domain_max);

            /// PoissonEq::PoissonVec<dim> poissonVec(octDA, &tree, DOF);
            /// poissonVec.setProblemDimensions(domain_min,domain_max);

            /// double * ux=octDA->getVecPointerToDof(uSolVecPtr,VAR::M_UI_U, false,false);
            /// double * frhs=octDA->getVecPointerToDof(uSolVecPtr,VAR::M_UI_F, false,false);
            /// double * Mfrhs=octDA->getVecPointerToDof(uSolVecPtr,VAR::M_UI_MF, false,false);
            double *ux = uSolVecPtr;
            /// double *frhs = fVecPtr;
            /// double *Mfrhs = mfVecPtr;
            double *dummy = dummyVecPtr;

            //TODO check if these interface are ready for 4D coordinates.

            /// std::function<void(double,double,double,double*)> f_rhs =[](const double x,const double y,const double z,double* var){
            ///     var[0]=1;
            /// };

            std::function<void(const double *, double*)> f_init =[](const double * xyz,double *var){
                var[0]=1;
            };

            octDA->setVectorByFunction(ux,f_init,false,false,DOF);
            /// octDA->setVectorByFunction(Mfrhs,f_init,false,false,DOF);
            /// octDA->setVectorByFunction(frhs,f_rhs,false,false,DOF);
            octDA->setVectorByFunction(dummy,f_init,false,false,DOF);

            /// poissonVec.computeVec(frhs,Mfrhs,1.0);

            // Warmup for matvec.
            for (int ii = 0; ii < numWarmup; ii++)
            {
                myPoissonMat.matVec(ux, dummy, 1.0);
            }

            // Clear the side effect of warmup.
            t_ghostexchange.clear();
            t_topdown.clear();
            t_bottomup.clear();
            t_treeinterior.clear();
            t_elemental.clear();
            t_matvec.clear();
            //resetAllTimers();
            // Benchmark the matvec.
            for (int ii = 0; ii < numRuns; ii++)
            {
                myPoissonMat.matVec(ux, dummy, 1.0);
            }

            /// double tol=1e-6;
            /// unsigned int max_iter=1000;
            /// myPoissonMat.cgSolve(ux,Mfrhs,max_iter,tol,0);

            /// const char * vNames[]={"m_uiU","m_uiFrhs","m_uiMFrhs"};
            /// octDA->vecTopvtu(uSolVecPtr,"poissonEq",(char**)vNames,false,false,DOF);
            /// octDA->destroyVector(uSolVec);

            delete octDA;
        }

    }


    void dump_profile_info(std::ostream& fout, const char *msgPrefix, double *params, const char **paramNames, unsigned int numParams, profiler_t* timers, const char ** names, unsigned int n ,MPI_Comm comm)
    {

        double stat;
        double stat_g[3*n];

        int rank, npes;
        MPI_Comm_rank(comm,&rank);
        MPI_Comm_size(comm,&npes);


        for(unsigned int i=0; i<n; i++)
        {
           stat=(timers[i].seconds); /// timers[i].num_calls ;     
           
           par::Mpi_Reduce(&stat,stat_g + 3*i + 0 ,1, MPI_MIN,0,comm);
           par::Mpi_Reduce(&stat,stat_g + 3*i + 1 ,1, MPI_SUM,0,comm);
           par::Mpi_Reduce(&stat,stat_g + 3*i + 2 ,1, MPI_MAX,0,comm);

           stat_g[ 3*i + 1] = stat_g[ 3*i + 1]/(double)npes;
 
        }

        par::Mpi_Reduce(&gRptSz.b1_treeSortSz,         &gDistRptSz.b1_treeSortSz,         1, MPI_SUM, 0, comm);
        par::Mpi_Reduce(&gRptSz.b1_treeConstructionSz, &gDistRptSz.b1_treeConstructionSz, 1, MPI_SUM, 0, comm);
        par::Mpi_Reduce(&gRptSz.b1_treeBalancingSz,    &gDistRptSz.b1_treeBalancingSz,    1, MPI_SUM, 0, comm);
        par::Mpi_Reduce(&gRptSz.b2_treeMatvecSz,       &gDistRptSz.b2_treeMatvecSz,       1, MPI_SUM, 0, comm);


        if(!rank)
        {
            fout << "msgPrefix\t" << "npes\t";
            for (unsigned int i = 0; i < numParams; i++)
            {
              fout << paramNames[i] << "\t";
            }
            fout << "treeSortSz\t" << "treeConstructionSz\t" << "treeBalancingSz\t" << "constrNumNodes\t"
                 << "treeMatvecSz\t" << "matvecNumNodes\t";
            for(unsigned int i=0; i<n; i++)
            {
               fout<<names[i]<<"(min)\t"<<names[i]<<"(mean)\t"<<names[i]<<"(max)\t";
            }

        }

        if(!rank)
            fout<<std::endl;

        if(!rank)
        {
            fout << msgPrefix << "\t";
            fout << npes << "\t";
            for (unsigned int i = 0; i < numParams; i++)
            {
              fout << params[i] << "\t";
            }
            fout << gDistRptSz.b1_treeSortSz << "\t"
                 << gDistRptSz.b1_treeConstructionSz << "\t"
                 << gDistRptSz.b1_treeBalancingSz << "\t"
                 << gDistRptSz.b1_globNodeSz << "\t"
                 << gDistRptSz.b2_treeMatvecSz << "\t"
                 << gDistRptSz.b2_globNodeSz << "\t";
            for(unsigned int i=0; i<n; i++)
            {
               fout<<stat_g[3*i + 0]<<"\t"<<stat_g[3*i + 1]<<"\t"<<stat_g[3*i+2]<<"\t";
            }
        }

        if(!rank)
            fout<<std::endl;

    }


}// end of namespace of bench



int main(int argc, char** argv)
{
    MPI_Init(&argc,&argv);
    DendroScopeBegin();

    signal(SIGABRT, handler);  // Register handler().

    int rank,npes;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm,&rank);
    MPI_Comm_size(comm,&npes);

    const unsigned int msgPrefixLimit = 64;
    
    if(argc<=1)
    {
        if(!rank)
            std::cout<<"usage :  "<<argv[0]<<" pts_per_core(weak scaling) maxdepth elementalOrder [lenPower2 (e.g. 4 for 16:1:1, default)] msgPrefix(<" << msgPrefixLimit << ")"<<std::endl;
        
        MPI_Abort(comm,0);
    }

    constexpr unsigned int dim = 3;

    const unsigned int pts_per_core = atoi(argv[1]);
    m_uiMaxDepth = atoi(argv[2]);
    const unsigned int eleOrder = atoi(argv[3]);

    char msgPrefix[2*msgPrefixLimit + 1];
    msgPrefix[0] = '\0';
    msgPrefix[msgPrefixLimit] = '\0';
    if (argc > 4 && argc <= 5)
      std::strncpy(msgPrefix, argv[4], msgPrefixLimit);
    else if (argc > 5)
      std::strncpy(msgPrefix, argv[5], msgPrefixLimit);

    int lengthPower2 = 4;
    if (argc > 5)
      lengthPower2 = atoi(argv[4]);

    _InitializeHcurve(dim);

    const unsigned int numTreeRuns = 1;
    const unsigned int numWarmup = 5;
    const unsigned int numRuns = 100;
    bench::bench_kernel<dim>(pts_per_core, numTreeRuns, numWarmup, numRuns, eleOrder, lengthPower2,comm);


    const char * param_names[] = {
        "pts_per_core",
        "eleOrder",
    };

    double params[] = {
        (double) pts_per_core,
        (double) eleOrder,
    };


    const char * counter_names[] = {
        "sort", 
        "constr", 
        "bal", 
        "adaptive_oda", 

        "matvec",
        "ghostexchange", 
        "topdown", 
        "bottomup", 
        /// "treeinterior", 
        "elemental", 
    };

    profiler_t counters[] = {
        bench::t_adaptive_tsort, 
        bench::t_adaptive_tconstr, 
        bench::t_adaptive_tbal, 
        bench::t_adaptive_oda, 

        bench::t_matvec,
        bench::t_ghostexchange, 
        bench::t_topdown, 
        bench::t_bottomup, 
        /// bench::t_treeinterior, 
        bench::t_elemental, 
    };

    bench::dump_profile_info(std::cout, msgPrefix, params,param_names,2, counters,counter_names,9, comm);

    _DestroyHcurve();
    DendroScopeEnd();
    MPI_Finalize();
    return 0;
}
