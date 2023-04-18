#include "hcurvedata.h"
#include "treeNode.h"
#include "oda.h"

#include <petsc.h>
#include <array>
#include <vector>

#include "sfcTreeLoop_matvec_io.h"

int main(int argc, char * argv[])
{
  constexpr int DIM = 2;

  PetscInitialize(&argc, &argv, NULL, NULL);
  DendroScopeBegin();
  _InitializeHcurve(DIM);

  int eleOrder = 2;
  int ndof = 1;
  m_uiMaxDepth = 10;
  /// int level = 2;
  int level = 5;

  ot::DA<DIM> octDA;
  std::array<unsigned int, DIM> a = {3,3};
  std::vector<ot::TreeNode<unsigned int, DIM>> treePart;
  ot::constructRegularSubdomainDA<DIM>(octDA, treePart, level, a, eleOrder, MPI_COMM_WORLD);

  std::vector<size_t> bdyIndex;
  octDA.getBoundaryNodeIndices(bdyIndex);

  std::cout << "octDA local num elements == " << octDA.getLocalElementSz() << "\n";
  std::cout << "octDA local num nodes == " << octDA.getLocalNodalSz() << "\n";
  std::cout << "octDA local bdyIndex.size() == " << bdyIndex.size() << "\n";


  // -------------------
  // Boundary Elements
  // -------------------

  /// std::vector<ot::TreeNode<unsigned int, DIM>> boundaryElementCoords;

  /// const bool visitEmpty = false;
  /// const unsigned int padlevel = 0;
  /// ot::MatvecBaseCoords<DIM> treeloop_coords(octDA.getTotalNodalSz(),
  ///                                           octDA.getElementOrder(),
  ///                                           visitEmpty,
  ///                                           padlevel,
  ///                                           octDA.getTNCoords(),
  ///                                           &(*treePart.cbegin()),
  ///                                           treePart.size(),
  ///                                           treePart.front(),
  ///                                           treePart.back());

  /// while (!treeloop_coords.isFinished())
  /// {
  ///   if (treeloop_coords.isPre() && treeloop_coords.isLeaf())
  ///   {
  ///     if (treeloop_coords.subtreeInfo().isElementBoundary())
  ///       boundaryElementCoords.push_back(treeloop_coords.getCurrentSubtree());

  ///     treeloop_coords.next();
  ///   }
  ///   else
  ///     treeloop_coords.step();
  /// }

  /// ot::printNodeCoords(&(*boundaryElementCoords.begin()),
  ///                     &(*boundaryElementCoords.cend()),
  ///                     eleOrder,
  ///                     std::cout);

  DendroScopeEnd();
  PetscFinalize();

  return 0;
}
