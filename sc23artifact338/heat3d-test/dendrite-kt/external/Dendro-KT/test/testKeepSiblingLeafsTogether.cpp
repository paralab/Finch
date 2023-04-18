
/**
 * @file testKeepSiblingLeafsTogether.cpp
 * @author Masado Ishii
 *
 * Test the utility function keepSiblingLeafsTogether().
 *
 * We will construct a complete tree based on a random Gaussian point cloud,
 * count the number of ranks with broken siblings. If there are children of
 * a sibling closer to the middle of the local partition, it doesn't count.
 */


#include "treeNode.h"
#include "tsort.h"
#include "parUtils.h"
#include "octUtils.h"
#include "hcurvedata.h"

#include <stdio.h>
#include <iostream>
#include <random>
#include <vector>



template <unsigned int dim>
bool testRandTree(MPI_Comm comm);


// ==============================
// main()
// ==============================
int main(int argc, char * argv[])
{
  MPI_Init(&argc, &argv);
  int totalSuccess = true;
  DendroScopeBegin();

  MPI_Comm comm = MPI_COMM_WORLD;

  int rProc, nProc;
  MPI_Comm_rank(comm, &rProc);
  MPI_Comm_size(comm, &nProc);

  const char * usageString = "Usage: %s dim\n";
  unsigned int inDim;

  if (argc < 2)
  {
    if (!rProc)
      printf(usageString, argv[0]);
    exit(1);
  }
  else
  {
    inDim   = strtol(argv[1], NULL, 0);
  }

  _InitializeHcurve(inDim);

  if (!rProc)
    printf("Test results: ");

  const char * resultColor;
  const char * resultName;

  // testRandTree
  int result_testRandTree, globResult_testRandTree;
  switch (inDim)
  {
    case 2: result_testRandTree = testRandTree<2>(comm); break;
    case 3: result_testRandTree = testRandTree<3>(comm); break;
    case 4: result_testRandTree = testRandTree<4>(comm); break;
    default: if (!rProc) printf("Dimension not supported.\n"); exit(1); break;
  }
  par::Mpi_Reduce(&result_testRandTree, &globResult_testRandTree, 1, MPI_SUM, 0, comm);
  totalSuccess = totalSuccess && !globResult_testRandTree;
  resultColor = globResult_testRandTree ? RED : GRN;
  resultName = globResult_testRandTree ? "FAILURE" : "success";
  if (!rProc)
    printf("\t[testRandTree](%s%s %d%s)", resultColor, resultName, globResult_testRandTree, NRM);

  if(!rProc)
    printf("\n");

  _DestroyHcurve();

  DendroScopeEnd();
  MPI_Finalize();

  return (!totalSuccess);
}





template <unsigned int dim>
bool testRandTree(MPI_Comm comm)
{
  // Test:
  // while mpirun -np <NP> ./tstKeepSiblingLeafsTogether <dim> > /dev/null ; do echo ; done

  using C = unsigned int;

  int rProc, nProc;
  MPI_Comm_rank(comm, &rProc);
  MPI_Comm_size(comm, &nProc);

  const unsigned int totalNumPoints = 50;
  const unsigned int numMyPoints = (totalNumPoints / nProc)
                                   + (rProc < totalNumPoints % nProc);
  const unsigned int numPrevPoints = (totalNumPoints / nProc) * rProc
                                     + (rProc < totalNumPoints % nProc ? rProc : totalNumPoints % nProc);

  unsigned int seed;
  int countInitial_glob = 0;
  int countFinal_glob = 0;
  std::vector<ot::TreeNode<C,dim>> tree;


  // Repeat until we get a test case where siblings are actually split.
  int trialSeeds = 0;
  const bool justOnce = false;
  do
  {
    trialSeeds++;
    tree.clear();

    if (!rProc)
    {
      // Pseudo-random number generators for point coordinates.
      std::random_device rd;
      seed = rd();
      // Can also set seed manually if needed.

      /// std::cerr << "Seed: " << seed << "\n";  // Wait till we know it's a good test.
    }
    MPI_Bcast(&seed, 1, MPI_UNSIGNED, 0, comm);
    std::mt19937 gen(seed);
    gen.discard((dim + 1) * numPrevPoints);

    std::uniform_int_distribution<C> coordDis(0, (1u << m_uiMaxDepth) - 1);
    std::uniform_int_distribution<unsigned int> levDis((m_uiMaxDepth + 1)/4, (m_uiMaxDepth + 1)/2);

    std::vector<ot::TreeNode<C,dim>> pointCoords;

    for (int ii = 0; ii < numMyPoints; ii++)
    {
      const unsigned int lev = levDis(gen);
      const C mask = (1u << (m_uiMaxDepth + 1)) - (1u << (m_uiMaxDepth - lev));

      std::array<C,dim> coords;
      for (int d = 0; d < dim; d++)
        coords[d] = coordDis(gen) & mask;

      pointCoords.emplace_back(1, coords, lev);
    }

    // Make points at least locally unique.
    ot::SFC_Tree<C,dim>::locTreeSort(
        &(*pointCoords.begin()),
        0,
        (ot::RankI) pointCoords.size(),
        1,
        m_uiMaxDepth,
        0);
    ot::SFC_Tree<C,dim>::locRemoveDuplicates(pointCoords);

    // Distributed tree construction.
    ot::SFC_Tree<C,dim>::distTreeConstruction(pointCoords, tree, 1, 0.0, comm);

    // Count number of procs with divided siblings, before filtering.
    {
      MPI_Comm nonemptys;
      MPI_Comm_split(comm, (tree.size() > 0 ? 1 : MPI_UNDEFINED), rProc, &nonemptys);
      if (tree.size())
      {
        countInitial_glob = ot::checkSiblingLeafsTogether(tree, nonemptys);
        MPI_Comm_free(&nonemptys);
      }
    }
  }
  while (!justOnce && trialSeeds < 1000 && !countInitial_glob);

  if (!rProc)
    std::cerr << "Seed: " << seed << " (trial " << trialSeeds << ")\n";  // Now we'll use the test.

  if (!rProc)
    std::cout << "countInitial_glob==" << countInitial_glob << " \n";

  ot::keepSiblingLeafsTogether<C,dim>(tree, comm);

  /// fprintf(stderr, "%*s[g%d] Finished keepSiblingLeafsTogether()\n", 40*rProc, "\0", rProc);

  /// //DEBUG
  /// for (int ii = 0; ii < tree.size(); ii++)
  /// {
  ///   fprintf(stdout, "%*c[%03d] == (%-17s).%u\n",
  ///       rProc * 40, ' ',
  ///       ii, tree[ii].getBase32Hex().data(), tree[ii].getLevel());
  /// }

  /// if (!tree.size())
  ///   fprintf(stdout, "%*s[g%d] I am empty!\n", 40*rProc, "\0", rProc);

  // Count number of procs with divided siblings, after filtering.
  {
    MPI_Comm nonemptys;
    MPI_Comm_split(comm, (tree.size() > 0 ? 1 : MPI_UNDEFINED), rProc, &nonemptys);
    if (tree.size())
    {
      countFinal_glob = ot::checkSiblingLeafsTogether(tree, nonemptys);
      MPI_Comm_free(&nonemptys);
    }
  }

  if (!rProc)
  {
    std::cout << "countFinal_glob==" << countFinal_glob << " \n";
  }

  return (countFinal_glob > 0);
}

