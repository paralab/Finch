/*
 * testTreeBalancing.cpp
 *   Test creation of balancing nodes and subsequent complete balanced tree.
 *
 * Masado Ishii  --  UofU SoC, 2019-01-10
 */


#include "treeNode.h"
#include "tsort.h"
#include "octUtils.h"

#include "hcurvedata.h"

#include "octUtils.h"
#include "oct2vtk.h"
#include <vector>

#include <assert.h>
#include <mpi.h>
#include <stdio.h>


// ...........................................................................
// ...........................................................................

//-------------------------------
// test_distOutputTreeBalancing()
//-------------------------------
void test_distOutputTreeBalancing(int numPoints, MPI_Comm comm = MPI_COMM_WORLD)
{
  int nProc, rProc;
  MPI_Comm_size(comm, &nProc);
  MPI_Comm_rank(comm, &rProc);

  using T = unsigned int;
  const unsigned int dim = 4;
  const unsigned int numChildren = 1u << dim;
  using TreeNode = ot::TreeNode<T,dim>;

  _InitializeHcurve(dim);

  std::vector<TreeNode> points = ot::getPts<T,dim>(numPoints);
  std::vector<TreeNode> treePart;

  const unsigned int maxPtsPerRegion = 32;

  const double loadFlexibility = 0.2;

  const T leafLevel = m_uiMaxDepth;
  const T firstVariableLevel = 1;      // Not sure about this whole root thing...

  printf("[%d] Starting distTreeBalancing()...\n", rProc);
  ot::SFC_Tree<T,dim>::distTreeBalancing(points, treePart, maxPtsPerRegion, loadFlexibility, comm);
  printf("[%d] Finished distTreeBalancing().\n", rProc);

  // Since dim==4 and we can only output octrees of dim 3 (or less?) to vtu,
  // use the slicing operator.
  unsigned int treeDepth = 0;
  for (const TreeNode &tn : treePart)
  {
    treeDepth = (tn.getLevel() > treeDepth ? tn.getLevel() : treeDepth);
  }
  unsigned int timeStep = 1u << (m_uiMaxDepth - treeDepth);
  const char dimNames[] = "XYZT";
  std::cout << "Output files will be placed in the directory '_output'.\n";
  for (unsigned int d = 0; d < dim; d++)
  {
    // Many slices, like a time series.
    for (unsigned int tIdx = 0, t = 0; t < (1u << m_uiMaxDepth); tIdx++, t += timeStep)
    {
      constexpr bool RM_DUPS_AND_ANC = false;
      constexpr bool RM_DUPS_ONLY = true;

      std::vector<ot::TreeNode<T,3>> slice3D;
      projectSliceKTree(&(*treePart.begin()), slice3D, (unsigned int) treePart.size(), d, (T) t);
      ot::SFC_Tree<T,3>::distRemoveDuplicates(slice3D, loadFlexibility, RM_DUPS_AND_ANC, comm);
      // Note that the partitioning of the slice is not related to partition of original tree.

      // Output to file with oct2vtu().
      char fPrefix[] =  "                                           ";  // beware buffer overflow.
      sprintf(fPrefix,  "_output/testSlice-%c-t%u", dimNames[d], tIdx);
      io::vtk::oct2vtu(&(*slice3D.begin()), (unsigned int) slice3D.size(), fPrefix, comm);
    }
  }
}



int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);
  DendroScopeBegin();

  int ptsPerProc = 200;
  if (argc > 1)
    ptsPerProc = strtol(argv[1], NULL, 0);

  test_distOutputTreeBalancing(ptsPerProc, MPI_COMM_WORLD);

  DendroScopeEnd();
  MPI_Finalize();

  return 0;
}
