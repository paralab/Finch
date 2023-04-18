
#include "treeNode.h"
#include "octUtils.h"
#include "parUtils.h"
#include "hcurvedata.h"

#include <stdio.h>
#include <iostream>
#include <vector>



template <unsigned int dim>
bool test(unsigned int lev, MPI_Comm comm);


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

  const char * usageString = "Usage: %s dim lev\n";
  unsigned int inDim;
  unsigned int inLev;

  if (argc < 3)
  {
    if (!rProc)
      printf(usageString, argv[0]);
    exit(1);
  }
  else
  {
    inDim   = strtol(argv[1], NULL, 0);
    inLev   = strtol(argv[2], NULL, 0);
  }

  _InitializeHcurve(inDim);

  const char * resultColor;
  const char * resultName;

  // test
  int result_test, globResult_test;
  switch (inDim)
  {
    case 2: result_test = test<2>(inLev, comm); break;
    case 3: result_test = test<3>(inLev, comm); break;
    case 4: result_test = test<4>(inLev, comm); break;
    default: if (!rProc) printf("Dimension not supported.\n"); exit(1); break;
  }
  par::Mpi_Reduce(&result_test, &globResult_test, 1, MPI_SUM, 0, comm);
  totalSuccess = totalSuccess && !globResult_test;
  resultColor = globResult_test ? RED : GRN;
  resultName = globResult_test ? "FAILURE" : "success";
  if (!rProc)
    printf("\t[test](%s%s %d%s)", resultColor, resultName, globResult_test, NRM);

  if(!rProc)
    printf("\n");

  _DestroyHcurve();

  DendroScopeEnd();
  MPI_Finalize();

  return (!totalSuccess);
}



template <unsigned int dim>
bool test(unsigned int lev, MPI_Comm comm)
{
  // Test:
  // while mpirun -np <NP> ./tstKeepSiblingLeafsTogether <dim> > /dev/null ; do echo ; done

  using C = unsigned int;

  int rProc, nProc;
  MPI_Comm_rank(comm, &rProc);
  MPI_Comm_size(comm, &nProc);

  std::vector<ot::TreeNode<C, dim>> treePart;

  // Create the tree.
  ot::createRegularOctree(treePart, lev, comm);

  // Check indeed all TreeNodes have level `lev' and volume is right.
  bool correctLevel = true;
  ot::RankI volume_loc = 0;
  ot::RankI volume_glob = 0;
  for (const auto &tn : treePart)
  {
    /// fprintf(stderr, "%*s[g%d] (%u)%-20s\n", 40*rProc, "\0", rProc,
    ///     tn.getLevel(), tn.getBase32Hex().data());

    if (tn.getLevel() == lev)
      volume_loc++;
    else
    {
      correctLevel = false;
      break;
    }
  }

  par::Mpi_Reduce(&volume_loc, &volume_glob, 1, MPI_SUM, 0, comm);

  /// fprintf(stderr, "%*s[g%d] correctLevel=%d, volume_loc==%llu\n",
  ///     40*rProc, "\0", rProc, 
  ///     int(correctLevel), volume_loc);

  return !(correctLevel && (rProc || volume_glob == (1u << (dim * lev))));
}

