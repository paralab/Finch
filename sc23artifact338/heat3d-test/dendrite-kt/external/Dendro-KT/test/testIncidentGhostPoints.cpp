#ifdef BUILD_WITH_PETSC
#include "petsc.h"
#endif

#include <iostream>

#include <treeNode.h>
#include <oda.h>
#include <sfcTreeLoop_matvec_io.h>
#include <filterFunction.h>

#include "genChannelPoints.h"

constexpr unsigned int DIM = 3;

/// #define CONDCOLOR(condition, color1, color2, body) \
///   ((condition) ? ( color1 body NRM ) : ( color2 body NRM ))

// ====================================
template <typename T>
class TypeWrapper
{
  private:
    T m_member;
  public:
    TypeWrapper(const T &t) : m_member(t) {}
    TypeWrapper(const TypeWrapper &other) = default;
    const T &get() const { return m_member; }
};


class NumLocalPoints : public TypeWrapper<size_t>
{
  using TypeWrapper<size_t>::TypeWrapper;
};
class NumGhostPoints : public TypeWrapper<size_t>
{
  using TypeWrapper<size_t>::TypeWrapper;
};
class NumIncidentGhostPoints : public TypeWrapper<size_t>
{
  using TypeWrapper<size_t>::TypeWrapper;
};


struct GhostIncidenceRate
{
  NumLocalPoints         m_numLocalPoints;
  NumGhostPoints         m_numGhostPoints;
  NumIncidentGhostPoints m_numIncidentGhostPoints;

  GhostIncidenceRate(const NumLocalPoints &numLocalPoints,
                     const NumGhostPoints &numGhostPoints,
                     const NumIncidentGhostPoints &numIncidentGhostPoints)
    : m_numLocalPoints(numLocalPoints),
      m_numGhostPoints(numGhostPoints),
      m_numIncidentGhostPoints(numIncidentGhostPoints)
  {}

  GhostIncidenceRate(const GhostIncidenceRate &other) = default;
  GhostIncidenceRate & operator=(const GhostIncidenceRate &other) = default;

  const NumLocalPoints & getNumLocalPoints() const { return m_numLocalPoints; }
  const NumGhostPoints & getNumGhostPoints() const { return m_numGhostPoints; }
  const NumIncidentGhostPoints & getNumIncidentGhostPoints() const { return m_numIncidentGhostPoints; }
};
// ====================================


// =====================================================================================================
GhostIncidenceRate countIncidentGhosts(const std::vector<ot::TreeNode<unsigned int, DIM>> &localOctList,
                                       const ot::DA<DIM> &da);
// =====================================================================================================





int main(int argc, char * argv[])
{
  // -----------------------------------------
#ifdef BUILD_WITH_PETSC
  PetscInitialize(&argc, &argv, NULL, NULL);
#else
  MPI_Init(&argc, &argv);
#endif
  DendroScopeBegin();
  _InitializeHcurve(DIM);
  // -----------------------------------------

  int npes, rank;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &rank);

  if (argc <= 1)
  {
    if(!rank)
    {
      std::cout << "usage :  " << argv[0] << " pts_per_core(weak scaling) maxdepth elementalOrder sfc_tol [lenPower2 (e.g. 4 for 16:1:1, default)]\n";
      std::flush(std::cout);
    }
    MPI_Abort(comm,0);
  }

  using T = unsigned int;
  using TreeNode = ot::TreeNode<T,DIM>;
  using OctList = std::vector<TreeNode>;
  using SFC_Tree = ot::SFC_Tree<T,DIM>;

  const unsigned int pts_per_core = atoi(argv[1]);
  m_uiMaxDepth = atoi(argv[2]);
  const unsigned int eleOrder = atoi(argv[3]);

  const double loadFlexibility = atof(argv[4]);

  int lengthPower2 = 4;
  if (argc > 5)
    lengthPower2 = atoi(argv[5]);

  const ibm::DomainDecider boxDecider = bench::getBoxDecider<DIM>(lengthPower2);

  std::vector<ot::TreeNode<unsigned int, DIM>> treePart = bench::getChannelPoints<DIM>(
      pts_per_core, lengthPower2, comm);
  SFC_Tree::distTreeSort(treePart, loadFlexibility, comm);

  ot::DistTree<unsigned int, DIM> distTree(treePart, comm);
  distTree.filterTree(boxDecider);

  ot::DA<DIM> da(distTree, comm, eleOrder, pts_per_core, loadFlexibility);

  const GhostIncidenceRate ghostIncidenceRate(
      countIncidentGhosts(distTree.getTreePartFiltered(), da)
      );
  const size_t numLocalNodes = ghostIncidenceRate.getNumLocalPoints().get();
  const size_t numGhostNodes = ghostIncidenceRate.getNumGhostPoints().get();
  const size_t numIncidentGhostNodes = ghostIncidenceRate.getNumIncidentGhostPoints().get();

  fprintf(stdout, "%03d\tsendTo\t%3d\trecvFrom\t%3d\t"
                  "numLocal==%lu\t"
                  "%snumGhost==%lu\t" NRM
                  "%snumIncidentGhost==%lu\t" NRM
                  "\n",
      rank,
      da.getNumOutboundRanks(), da.getNumInboundRanks(),
      numLocalNodes,
      (numGhostNodes < 0.5 * numLocalNodes ? GRN : RED ), numGhostNodes,
      (numIncidentGhostNodes >= 0.9 * numGhostNodes ? GRN : RED ), numIncidentGhostNodes
      );

  // -----------------------------------------
  _DestroyHcurve();
  DendroScopeEnd();
#ifdef BUILD_WITH_PETSC
  PetscFinalize();
#else
  MPI_Finalize();
#endif
  // -----------------------------------------

  return 0;
}


GhostIncidenceRate countIncidentGhosts(const std::vector<ot::TreeNode<unsigned int, DIM>> &localOctList,
                                       const ot::DA<DIM> &da)
{
  const size_t totalVectorSz = da.getTotalNodalSz();
  const size_t localVectorSz = da.getLocalNodalSz();
  const size_t numGhostNodes = totalVectorSz - localVectorSz;

  std::vector<int> incidence(totalVectorSz, 0);

  const unsigned int ndofs = 1;
  const bool visitEmpty = false;
  const unsigned int padLevel = 0;

  ot::MatvecBaseOut<DIM, int, true> mvecBaseOut(
      totalVectorSz,
      ndofs,
      da.getElementOrder(),
      visitEmpty,
      padLevel,
      da.getTNCoords(),
      &(*localOctList.cbegin()),
      localOctList.size(),
      localOctList.front(),
      localOctList.back()
      );

  const std::vector<int> elementNodesAreOne(da.getNumNodesPerElement(), 1);

  while (!mvecBaseOut.isFinished())
  {
    if (mvecBaseOut.isLeaf() && !mvecBaseOut.isPre())
      mvecBaseOut.subtreeInfo().overwriteNodeValsOut(&(*elementNodesAreOne.cbegin()));
    mvecBaseOut.step();

    // Note it is also possible to write during Pre&&Leaf
    // and then use next() to skip the Post&&Leaf phase.
  }
  mvecBaseOut.finalize(&(*incidence.begin()));

  size_t numIncidentGhostPoints = 0;
  for (size_t ii = 0; ii < da.getLocalNodeBegin(); ++ii)
    if (incidence[ii] > 0)
      numIncidentGhostPoints++;
  for (size_t ii = da.getLocalNodeBegin() + localVectorSz; ii < totalVectorSz; ++ii)
    if (incidence[ii] > 0)
      numIncidentGhostPoints++;

  return GhostIncidenceRate(NumLocalPoints(localVectorSz),
                            NumGhostPoints(numGhostNodes),
                            NumIncidentGhostPoints(numIncidentGhostPoints));
}


