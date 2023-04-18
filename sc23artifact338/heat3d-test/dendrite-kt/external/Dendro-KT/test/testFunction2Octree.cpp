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


template <typename C, unsigned int dim>
std::vector<ot::TreeNode<C,dim>> getTree(double wavelet_tol, double partition_tol, unsigned int eOrder, MPI_Comm comm)
{
  unsigned int varIndex[1] = {0};

  // Function for generating the tree, as in exHeatEq.
  double g_min = 0.0;
  double g_max = 1.0;
  double d_min = -0.5;
  double d_max =  0.5;
  double Rg = g_max - g_min;
  double Rd = d_max - d_min;
  std::function<void(const double *, double*)> f_rhs = [d_min, d_max, g_min, g_max, Rg, Rd](const double *x, double *var)
  {
    var[0] = -12*M_PI*M_PI;
    for (unsigned int d = 0; d < dim; d++)
      var[0] *= sin(2*M_PI*(((x[d]-g_min)/Rg)*Rd+d_min));
  };

  std::vector<ot::TreeNode<C,dim>> completeTree, balancedTree;
  ot::function2Octree<C,dim>(f_rhs, 1, varIndex, 1, completeTree, m_uiMaxDepth, wavelet_tol, partition_tol, eOrder, comm);

  balancedTree = completeTree;
  ot::SFC_Tree<C,dim>::distTreeBalancing(completeTree, balancedTree, 1, partition_tol, comm);

  return balancedTree;
}

template <typename C, unsigned int dim>
bool distCompareTrees(const std::vector<ot::TreeNode<C,dim>> &treeA, const std::vector<ot::TreeNode<C,dim>> &treeB, MPI_Comm comm)
{
  int equalLoc = true, equalGlob = false;
  if (treeA.size() != treeB.size())
    equalLoc = false;
  else
    for (unsigned int ii = 0; ii < treeA.size(); ii++)
      if (treeA[ii] != treeB[ii])
      {
        equalLoc = false;
        break;
      }

  par::Mpi_Allreduce(&equalLoc, &equalGlob, 1, MPI_SUM, comm);
  return equalGlob;
}


// ==============================================================
// main_(): Implementation after parsing, getting dimension, etc.
// ==============================================================
template <unsigned int dim>
int main_ (Parameters &pm, MPI_Comm comm)
{
  using C = unsigned int;
  const unsigned int m_uiDim = dim;

  int rank, npes;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &npes);

  m_uiMaxDepth = pm.maxDepth;
  const double wavelet_tol = pm.waveletTol;
  const double partition_tol = pm.partitionTol;
  const unsigned int eOrder = pm.eleOrder;

  // Create equivalence classes of results from calling getTree() (function2Octree()).
  std::vector<std::vector<ot::TreeNode<C, dim>>> treeList(npes);
  std::vector<std::vector<int>> subsets;
  for (int n = 1; n <= npes; n++)
  {
    if (!rank)
      printf("Starting with subcomm n==%d\n", n);

    MPI_Comm subcomm;
    par::splitCommUsingSplittingRank(n, &subcomm, comm);

    std::vector<ot::TreeNode<C, dim>> tree;
    if (rank < n)
      tree = getTree<C, dim>(wavelet_tol, partition_tol, eOrder, subcomm);

    ot::SFC_Tree<C,dim>::distTreeSort(tree, 0.0, comm);

    int matchingRepresentative = 0;   // 0 means no match.
    auto subsetsIter = subsets.begin();
    for ( ; subsetsIter < subsets.end(); subsetsIter++)
    {
      int rep = subsetsIter->front();
      if (distCompareTrees<C,dim>(tree, treeList[rep-1], comm))
      {
        matchingRepresentative = rep;
        break;
      }
    }

    if (matchingRepresentative == 0)
    {
      std::swap(treeList[n-1], tree);
      subsets.emplace_back(1, n);
    }
    else
    {
      subsetsIter->push_back(n);
    }
  }

  // Print equivalence classes.
  if (!rank)
  {
    printf("Equivalence classes (%s%lu%s)\n", (subsets.size() > 1 ? RED : GRN), subsets.size(), NRM);
    for (auto &&s : subsets)
    {
      printf("(");
      for (int x : s)
        printf(" %d ", x);
      printf(")");
      printf("\n");
    }
  }

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
  MPI_Init(&argc, &argv);
  int returnCode = 1;
  DendroScopeBegin();

  int rProc, nProc;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &rProc);
  MPI_Comm_size(comm, &nProc);

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

  DendroScopeEnd();
  MPI_Finalize();

  return returnCode;
}
