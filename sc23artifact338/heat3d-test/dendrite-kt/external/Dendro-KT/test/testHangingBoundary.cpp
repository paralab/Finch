#include "hcurvedata.h"
#include "treeNode.h"
#include "oda.h"

#include <petsc.h>
#include <array>
#include <vector>

#include "sfcTreeLoop_matvec_io.h"
#include "mathUtils.h"

int main(int argc, char * argv[])
{
  constexpr int DIM = 3;

  PetscInitialize(&argc, &argv, NULL, NULL);
  DendroScopeBegin();
  _InitializeHcurve(DIM);

  int nProc, rProc;
  MPI_Comm_size(MPI_COMM_WORLD, &nProc);
  MPI_Comm_rank(MPI_COMM_WORLD, &rProc);

  int eleOrder = 1;
  int ndof = 1;
  m_uiMaxDepth = 10;
  int level = 5;

  std::array<double, DIM> physDims = {1.0, 0.5, 0.5};

  std::vector<ot::TreeNode<unsigned int, DIM>> treePart;
  std::vector<ot::TreeNode<unsigned int, DIM>> treeCoarse;
  const ot::TreeNode<unsigned int, DIM> treeRoot;
  if (rProc == 0)
  {
    treeCoarse.push_back(treeRoot.getChildMorton(0));

    const ot::TreeNode<unsigned int, DIM> child1 = treeRoot.getChildMorton(1);
    treeCoarse.push_back(child1.getChildMorton(0));
    treeCoarse.push_back(child1.getChildMorton(1));
    treeCoarse.push_back(child1.getChildMorton(2));
    treeCoarse.push_back(child1.getChildMorton(3));
    treeCoarse.push_back(child1.getChildMorton(4));
    treeCoarse.push_back(child1.getChildMorton(5));
    treeCoarse.push_back(child1.getChildMorton(6));
    treeCoarse.push_back(child1.getChildMorton(7));

    const unsigned int numChildren = (1u << DIM);
    for (const ot::TreeNode<unsigned int, DIM> &tn : treeCoarse)
      for (int childM = 0; childM < numChildren; ++childM)
        treePart.push_back(tn.getChildMorton(childM));

    ot::SFC_Tree<unsigned int, DIM>::locTreeSort(treePart);
  }

  ot::DistTree<unsigned int, DIM> distTree(treePart, MPI_COMM_WORLD);
  distTree.filterTree((typename ot::DistTree<unsigned int, DIM>::BoxDecider)(physDims));
  treePart = distTree.getTreePartFiltered();
  ot::DA<DIM> octDA(distTree, MPI_COMM_WORLD, eleOrder, 1, 0.3);

  fprintf(stdout, "[rank%d] treePart.size() == %lu\n", rProc, treePart.size());
  fprintf(stdout, "[rank%d] octDA.getLocalNodalSz() == %lu\n", rProc, octDA.getLocalNodalSz());
  fprintf(stdout, "[rank%d] octDA.getGlobalNodeSz() == %llu\n", rProc, octDA.getGlobalNodeSz());

  /// for (size_t ii = 0; ii < octDA.getTotalNodalSz(); ++ii)
  ///   ot::printtn(octDA.getTNCoords()[ii], 2, std::cout) << "\n";

  std::vector<std::array<double, DIM>> bdryPoints;

  const bool visitEmpty = false;
  const unsigned int padlevel = 0;
  ot::MatvecBaseCoords<DIM> treeloop_coords(octDA.getTotalNodalSz(),
                                            octDA.getElementOrder(),
                                            visitEmpty,
                                            padlevel,
                                            octDA.getTNCoords(),
                                            &(*treePart.cbegin()),
                                            treePart.size(),
                                            treePart.front(),
                                            treePart.back());

  const size_t npe = intPow((1+eleOrder), DIM);
  while (!treeloop_coords.isFinished())
  {
    if (treeloop_coords.isPre() && treeloop_coords.isLeaf())
    {
      if (treeloop_coords.subtreeInfo().isElementBoundary())
      {
        const std::vector<bool> & nodeBoundary = treeloop_coords.subtreeInfo().getLeafNodeBdry();
        for (size_t nidx = 0; nidx < npe; ++nidx)
        {
          std::array<double, DIM> physCoords;
          for (int d = 0; d < DIM; ++d)
            physCoords[d] = treeloop_coords.subtreeInfo().getNodeCoords()[DIM * nidx + d];

          if (nodeBoundary[nidx])
            bdryPoints.push_back(physCoords);
        }
      }

      treeloop_coords.next();
    }
    else
      treeloop_coords.step();
  }

  std::sort(bdryPoints.begin(), bdryPoints.end());
  auto last = std::unique(bdryPoints.begin(), bdryPoints.end());
  bdryPoints.erase(last, bdryPoints.end());

  const size_t numTotalBdryNodes = bdryPoints.size();

  std::ofstream outFile("tstCoords.txt");

  size_t resetCount = 0;
  double resetX = 0.0;
  double resetY = 0.0;
  for (const auto &pt : bdryPoints)
  {
    if (pt[0] != resetX)
    {
      std::cout << "---------------\n  (" << resetCount << ")\n\n";
      resetCount = 0;
    }
    resetCount++;
    resetX = pt[0];

    std::cout << "[" << pt[0] << " " << pt[1] << " " << pt[2] << "]";
    outFile << pt[0] << " " << pt[1] << " " << pt[2] << "\n";

    if (pt[0] == 0.5
        && !((pt[1] == 0.0 || pt[1] == 0.5) && (pt[2] == 0.0 || pt[2] == 0.25 || pt[2] == 0.5)
            || (pt[2] == 0.0 || pt[2] == 0.5) && (pt[1] == 0.0 || pt[1] == 0.25 || pt[1] == 0.5) )
       )
      std::cout << " *";
    
    std::cout << "\n";
  }
  std::cout << "---------------\n  (" << resetCount << ")\n\n";

  outFile.close();

  fprintf(stdout, "[rank%d] numTotalBdryNodes == %lu\n", rProc, numTotalBdryNodes);

  DendroScopeEnd();
  PetscFinalize();

  return 0;
}
