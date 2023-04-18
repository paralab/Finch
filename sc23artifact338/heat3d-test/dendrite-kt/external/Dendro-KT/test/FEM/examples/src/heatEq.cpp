//
// Created by milinda on 11/21/18.
// Modified by masado on 04/24/19.
//


#include "treeNode.h"
#include "mpi.h"
#include "tsort.h"
#include "dendro.h"
#include "octUtils.h"
#include "functional"
#include "refel.h"
#include "heatMat.h"
#include "heatVec.h"

#ifdef BUILD_WITH_PETSC
  #include <petsc.h>
  #include <petscvec.h>
  #include <petscksp.h>
#endif

// =======================================================
// Parameters: Change these and the options in get_args().
// =======================================================
struct Parameters
{
  unsigned int dim;
  unsigned int maxDepth;
  double waveletTol;
  double partitionTol;
  unsigned int eleOrder;
};
// =======================================================


// ==============================================================
// main_(): Implementation after parsing, getting dimension, etc.
// ==============================================================
template <unsigned int dim>
int main_ (Parameters &pm, MPI_Comm comm)
{
    const unsigned int m_uiDim = dim;

    int rank, npes;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &npes);

    m_uiMaxDepth = pm.maxDepth;
    const double wavelet_tol = pm.waveletTol;
    const double partition_tol = pm.partitionTol;
    const unsigned int eOrder = pm.eleOrder;

    double tBegin = 0, tEnd = 10, th = 0.01;

    RefElement refEl(m_uiDim,eOrder);

    /// Point<dim> grid_min(0, 0, 0);
    /// Point<dim> grid_max(1, 1, 1);

    /// Point<dim> domain_min(-0.5,-0.5,-0.5);
    /// Point<dim> domain_max(0.5,0.5,0.5);

    /// double Rg_x=(grid_max.x()-grid_min.x());
    /// double Rg_y=(grid_max.y()-grid_min.y());
    /// double Rg_z=(grid_max.z()-grid_min.z());

    /// double Rd_x=(domain_max.x()-domain_min.x());
    /// double Rd_y=(domain_max.y()-domain_min.y());
    /// double Rd_z=(domain_max.z()-domain_min.z());

    /// const Point<dim> d_min=domain_min;
    /// const Point<dim> d_max=domain_max;

    /// const Point<dim> g_min=grid_min;
    /// const Point<dim> g_max=grid_max;

    // For now must be anisotropic.
    double g_min = 0.0;
    double g_max = 1.0;
    double d_min = -0.5;
    double d_max =  0.5;
    double Rg = g_max - g_min;
    double Rd = d_max - d_min;
    const Point<dim> domain_min(d_min, d_min, d_min);
    const Point<dim> domain_max(d_max, d_max, d_max);

    /// std::function<void(const double *, double*)> f_rhs =[d_min,d_max,g_min,g_max,Rg_x,Rg_y,Rg_z,Rd_x,Rd_y,Rd_z](const double *x, double* var){
    ///     var[0]=(-12*M_PI*M_PI*sin(2*M_PI*(((x[0]-g_min.x())/(Rg_x))*(Rd_x)+d_min.x()))*sin(2*M_PI*(((x[1]-g_min.y())/(Rg_y))*(Rd_y)+d_min.y()))*sin(2*M_PI*(((x[2]-g_min.z())/(Rg_z))*(Rd_z)+d_min.z())));
    ///     //var[1]=(-12*M_PI*M_PI*sin(2*M_PI*(((x[0]-g_min.x())/(Rg_x))*(Rd_x)+d_min.x()))*sin(2*M_PI*(((x[1]-g_min.y())/(Rg_y))*(Rd_y)+d_min.y()))*sin(2*M_PI*(((x[2]-g_min.z())/(Rg_z))*(Rd_z)+d_min.z())));
    ///     //var[2]=(-12*M_PI*M_PI*sin(2*M_PI*(((x[0]-g_min.x())/(Rg_x))*(Rd_x)+d_min.x()))*sin(2*M_PI*(((x[1]-g_min.y())/(Rg_y))*(Rd_y)+d_min.y()))*sin(2*M_PI*(((x[2]-g_min.z())/(Rg_z))*(Rd_z)+d_min.z())));
    /// };
    std::function<void(const double *, double*)> f_rhs = [d_min, d_max, g_min, g_max, Rg, Rd](const double *x, double *var)
    {
      var[0] = -12*M_PI*M_PI;
      for (unsigned int d = 0; d < dim; d++)
        var[0] *= sin(2*M_PI*(((x[d]-g_min)/Rg)*Rd+d_min));
    };

    std::function<void(const double *, double*)> f_init =[/*d_min,d_max,g_min,g_max,Rg_x,Rg_y,Rg_z,Rd_x,Rd_y,Rd_z*/](const double *x, double *var){
        var[0]=0;//(-12*M_PI*M_PI*sin(2*M_PI*(((x[0]-g_min.x())/(Rg_x))*(Rd_x)+d_min.x()))*sin(2*M_PI*(((x[1]-g_min.y())/(Rg_y))*(Rd_y)+d_min.y()))*sin(2*M_PI*(((x[2]-g_min.z())/(Rg_z))*(Rd_z)+d_min.z())));
        //var[1]=0;
        //var[2]=0;
    };


    ot::DistTree<unsigned, dim> distTree =
        ot::DistTree<unsigned, dim>::constructDistTreeByFunc(f_rhs, 1, comm, eOrder, wavelet_tol, partition_tol);
    ot::DA<dim> *octDA = new ot::DA<dim>(distTree, comm, eOrder, 100, partition_tol);
    const std::vector<ot::TreeNode<unsigned, dim>> &treePart = distTree.getTreePartFiltered();
    assert(treePart.size() > 0);

#ifndef BUILD_WITH_PETSC
    //
    // Non-PETSc version.
    //

    // There are three vectors that happen to have the same sizes but are logically separate.
    std::vector<double> ux, frhs, Mfrhs;
    octDA->createVector(ux, false, false, 1);
    octDA->createVector(frhs, false, false, 1);
    octDA->createVector(Mfrhs, false, false, 1);

    HeatEq::HeatMat<dim> heatMat(octDA, &treePart,1);
    heatMat.setProblemDimensions(domain_min,domain_max);

    HeatEq::HeatVec<dim> heatVec(octDA, &treePart,1);
    heatVec.setProblemDimensions(domain_min,domain_max);

    octDA->setVectorByFunction(ux.data(),f_init,false,false,1);
    octDA->setVectorByFunction(Mfrhs.data(),f_init,false,false,1);
    octDA->setVectorByFunction(frhs.data(),f_rhs,false,false,1);

    heatVec.computeVec(&(*frhs.cbegin()), &(*Mfrhs.begin()), 1.0);


    double tol=1e-6;
    unsigned int max_iter=1000;
    heatMat.cgSolve(&(*ux.begin()), &(*Mfrhs.begin()), max_iter, tol);

    // TODO
    // octDA->vecTopvtu(...);

    octDA->destroyVector(ux);
    octDA->destroyVector(frhs);
    octDA->destroyVector(Mfrhs);

#else
    //
    // PETSc version.
    //

    // There are three vectors that happen to have the same sizes but are logically separate.
    Vec ux, frhs, Mfrhs;
    octDA->petscCreateVector(ux, false, false, 1);
    octDA->petscCreateVector(frhs, false, false, 1);
    octDA->petscCreateVector(Mfrhs, false, false, 1);

    HeatEq::HeatMat<dim> heatMat(octDA, &treePart,1);
    heatMat.setProblemDimensions(domain_min,domain_max);

    HeatEq::HeatVec<dim> heatVec(octDA, &treePart,1);
    heatVec.setProblemDimensions(domain_min,domain_max);

    octDA->petscSetVectorByFunction(ux, f_init, false, false, 1);
    octDA->petscSetVectorByFunction(Mfrhs, f_init, false, false, 1);
    octDA->petscSetVectorByFunction(frhs, f_rhs, false, false, 1);

    heatVec.computeVec(frhs, Mfrhs, 1.0);

    double tol=1e-6;
    unsigned int max_iter=1000;

    Mat matrixFreeMat;
    heatMat.petscMatCreateShell(matrixFreeMat);

    // PETSc solver context.
    KSP ksp;
    PetscInt numIterations;

    KSPCreate(comm, &ksp);
    KSPSetOperators(ksp, matrixFreeMat, matrixFreeMat);
    KSPSetTolerances(ksp, tol, PETSC_DEFAULT, PETSC_DEFAULT, max_iter);

    KSPSolve(ksp, Mfrhs, ux);
    KSPGetIterationNumber(ksp, &numIterations);

    if (!rank)
      std::cout << " finished at iteration " << numIterations << " ...\n";

    KSPDestroy(&ksp);

    // Now that we have an approximate solution, test convergence by evaluating the residual.
    Vec residual;
    octDA->petscCreateVector(residual, false, false, 1);
    heatMat.matVec(ux, residual);
    VecAXPY(residual, -1.0, Mfrhs);
    PetscScalar normr, normb;
    VecNorm(Mfrhs, NORM_INFINITY, &normb);
    VecNorm(residual, NORM_INFINITY, &normr);
    PetscScalar rel_resid_err = normr / normb;

    if (!rank)
      std::cout << "Final relative residual error == " << rel_resid_err << ".\n";

    // TODO
    // octDA->vecTopvtu(...);

    octDA->petscDestroyVec(ux);
    octDA->petscDestroyVec(frhs);
    octDA->petscDestroyVec(Mfrhs);
    octDA->petscDestroyVec(residual);

#endif

    if(!rank)
        std::cout<<" end of heatEq: "<<std::endl;

    delete octDA;

    /// MPI_Finalize();
    return 0;
}
// ==============================================================


//
// get_args()
//
bool get_args(int argc, char * argv[], Parameters &pm, MPI_Comm comm)
{
  int rProc, nProc;
  MPI_Comm_rank(comm, &rProc);
  MPI_Comm_size(comm, &nProc);

  // ========================
  // Set up accepted options.
  // ========================
  enum CmdOptions                           { progName, opDim, opMaxDepth, opWaveletTol, opPartitionTol, opEleOrder, NUM_CMD_OPTIONS };
  const char *cmdOptions[NUM_CMD_OPTIONS] = { argv[0], "dim", "maxDepth", "waveletTol", "partitionTol", "eleOrder", };
  const unsigned int firstOptional = NUM_CMD_OPTIONS;  // All required.
  // ========================

  // Fill argVals.
  std::array<const char *, NUM_CMD_OPTIONS> argVals;
  argVals.fill("");
  for (unsigned int op = 0; op < argc; op++)
    argVals[op] = argv[op];

  // Check if we have the required arguments.
  if (argc < firstOptional)
  {
    if (!rProc)
    {
      std::cerr << "Usage: ";
      unsigned int op = 0;
      for (; op < firstOptional; op++)
        std::cerr << cmdOptions[op] << " ";
      for (; op < NUM_CMD_OPTIONS; op++)
        std::cerr << "[" << cmdOptions[op] << "] ";
      std::cerr << "\n";
    }
    return false;
  }

  // ================
  // Parse arguments.
  // ================
  pm.dim      = static_cast<unsigned int>(strtoul(argVals[opDim], NULL, 0));
  pm.maxDepth = static_cast<unsigned int>(strtoul(argVals[opMaxDepth], NULL, 0));
  pm.eleOrder = static_cast<unsigned int>(strtoul(argVals[opEleOrder], NULL, 0));
  pm.waveletTol   = strtod(argVals[opWaveletTol], NULL);
  pm.partitionTol = strtod(argVals[opPartitionTol], NULL);
  // ================

  // Replay arguments.
  constexpr bool replayArguments = true;
  if (replayArguments && !rProc)
  {
    for (unsigned int op = 1; op < NUM_CMD_OPTIONS; op++)
      std::cout << YLW << cmdOptions[op] << "==" << argVals[op] << NRM << " \n";
    std::cout << "\n";
  }

  return true;
}


//
// main()
//
int main(int argc, char * argv[])
{
#ifndef BUILD_WITH_PETSC
  MPI_Init(&argc, &argv);
#else
  PetscInitialize(&argc, &argv, NULL, NULL);
#endif

  int rProc, nProc;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &rProc);
  MPI_Comm_size(comm, &nProc);

  int returnCode = 1;

  Parameters pm;
  unsigned int &dim = pm.dim;
  if (get_args(argc, argv, pm, comm))
  {
    int synchronize;
    MPI_Bcast(&synchronize, 1, MPI_INT, 0, comm);

    _InitializeHcurve(dim);

    // Convert dimension argument to template parameter.
    switch(dim)
    {
      case 2: returnCode = main_<2>(pm, comm); break;
      case 3: returnCode = main_<3>(pm, comm); break;
      case 4: returnCode = main_<4>(pm, comm); break;
      default:
        if (!rProc)
          std::cerr << "Dimension " << dim << " not currently supported.\n";
    }

    _DestroyHcurve();
  }

#ifndef BUILD_WITH_PETSC
  MPI_Finalize();
#else
  PetscFinalize();
#endif

  return returnCode;
}
