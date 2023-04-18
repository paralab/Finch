#ifdef BUILD_WITH_PETSC
#include "petsc.h"
#endif
#include <iostream>
#include <distTree.h>
#include <oda.h>
#include <point.h>
#include <sfcTreeLoop_matvec_io.h>
#include <octUtils.h>

#include "oct2vtk.h"

ibm::Partition DomainDecider(const double *physCoords, double physSize);

constexpr unsigned int DIM = 3;
typedef ot::TreeNode<unsigned int, DIM> TREENODE;
constexpr unsigned int nchild = 1u << DIM;
static constexpr double RADIUS = 0.45;
static constexpr double scaling = 16.0; /// The octant Coords are scaled by 16 in each direction.
/// The three spheres in the domain
static double centerSphere[3][3]{{0.5, 0.5, 8.0},
                                 {0.5, 0.0, 8.0 + RADIUS * 2 + 0.01},
                                 {0.5, 1.0, 8.0 + RADIUS * 2 + 0.01}};
static const double DomainExt[3] = {1.0, 1.0, 16.0}; /// Physical domain of 1 X 1 X 16 is carved out from [16,16,16]
static unsigned int maxLevel;
enum CREATION_STAGE : bool {
    INITIAL = true,
    FINAL = false
};
static CREATION_STAGE creationStage = CREATION_STAGE::INITIAL;

int dbgRank = 0;

void checkHangingNodes(ot::DA<DIM>* octDA, const std::vector<TREENODE> & treePart){
  const size_t sz = octDA->getTotalNodalSz();
  auto partFront = octDA->getTreePartFront();
  auto partBack = octDA->getTreePartBack();
  const auto tnCoords = octDA->getTNCoords();
  double counter = 0;

  const std::vector<int> nodeExists(octDA->getTotalNodalSz(), 1);

  std::ofstream fout("Hanging" + std::to_string(octDA->getRankActive()) +".txt");
  std::vector<double> hangingElements(octDA->getLocalElementSz(),0.0);
  ot::MatvecBaseIn<DIM,int,false> loop(sz,1,octDA->getElementOrder(), false,0,tnCoords, &(*nodeExists.cbegin()),&(*treePart.cbegin()),treePart.size(),*partFront,*partBack);

  const int rankActive = octDA->getRankActive();

  bool atLeastOneWrong = false;

  while(!loop.isFinished()){
    if (loop.isPre() && loop.subtreeInfo().isLeaf()) {
      if(loop.subtreeInfo().isElementBoundary()){
        const int * existence = loop.subtreeInfo().readNodeValsIn();
        const double *coords = loop.subtreeInfo().getNodeCoords();

        bool hangingProblem = false;

        for(int i = 0; i < 8; i++){
          if(existence[i] == 0){ // This should never happen if all interpolation nodes exist.
            hangingProblem = true;
          }
        }

        if (hangingProblem)
        {
          hangingElements[counter] = 1.0;
          for(int j = 0; j < 8; j++) {
            /// fout << coords[DIM*j+0]*scaling << " " << coords[DIM*j+1]*scaling << " "<< coords[DIM*j+2]*scaling << "\n";
            fout << coords[DIM*j+0] << " " << coords[DIM*j+1] << " "<< coords[DIM*j+2] << "\n";
          }
          fout << "\n";
          printf("%*sSomething wrong happened, rank(%d). level==%d, child_m==%d, myVals==[%d %d %d %d %d %d %d %d]\n",
              5*rankActive, "", rankActive,
              loop.subtreeInfo().getCurrentSubtree().getLevel(),
              loop.subtreeInfo().getCurrentSubtree().getMortonIndex(),
              existence[0],
              existence[1],
              existence[2],
              existence[3],
              existence[4],
              existence[5],
              existence[6],
              existence[7]
              );
        }

      }
      counter++;
      loop.next();
    }
    else{
      loop.step();
    }
  }
  std::cout << counter << " "<< octDA->getLocalElementSz() << "\n";
  fout.close();
}
template<unsigned int dim>
void printTree(const std::vector<ot::TreeNode<unsigned int, dim>> &treePart, int elev) {
  std::cout << "Tree\n";
  for (const ot::TreeNode<unsigned int, dim> &tn : treePart) {
    ot::printtn(tn, elev, std::cout);
    std::cout << "\n";
  }
  std::cout << "\n";
}


/**
 * This function only refines the elements marked as boundary octants
 * @param octDA
 * @param refineFlags
 */

bool generateRefinementFlags(const std::vector<TREENODE> &treePart,
                             std::vector<ot::OCT_FLAGS::Refine> &refineFlags) {

  bool refineAny = false;

  refineFlags.resize(treePart.size());
  std::fill(refineFlags.begin(), refineFlags.end(), ot::OCT_FLAGS::Refine::OCT_NO_CHANGE);

  size_t reffIdx = 0;

  for (const TREENODE &tn : treePart)
  {
    std::array<double, nchild*DIM> physCoords;
    double physSize;
    ot::treeNode2Physical(tn, &physCoords[0*DIM+0], physSize);

    for (int child = 1; child < nchild; ++child)
      for (int axis = 0; axis < DIM; ++axis)
        physCoords[child*DIM + axis] = physCoords[0*DIM + axis] + bool(child & (1u << axis)) * physSize;

    for (double &c : physCoords)
      c *= scaling;

    /// Refining around the sphere:
    // Logic: If any one of the point of the element is inside the sphere Refine
    for (int nSphere = 0; nSphere < 3; nSphere++) {
      for (int numCoords = 0; numCoords < nchild; numCoords++) {
        double distX = centerSphere[nSphere][0] - physCoords[numCoords * DIM + 0];
        double distY = centerSphere[nSphere][1] - physCoords[numCoords * DIM + 1];
        double distZ = centerSphere[nSphere][2] - physCoords[numCoords * DIM + 2];
        double dist = distX * distX + distY * distY + distZ * distZ;
        if (dist < (RADIUS * RADIUS) and tn.getLevel() < maxLevel) {
          refineFlags[reffIdx] = ot::OCT_FLAGS::Refine::OCT_REFINE;
          refineAny = true;
          break;
        }
      }
    }

    reffIdx++;
  }

  return refineAny;
}
                             

ibm::Partition DomainDecider(const double *physCoords, double physSize)
{
  double coords[3];
  for (int ii = 0; ii < 3; ++ii)
    coords[ii] = physCoords[ii] * scaling;
  const double size = physSize * scaling;

  const bool isOutCube =
      (0.0 < coords[0] and coords[0] + size < DomainExt[0]) and
      (0.0 < coords[1] and coords[1] + size < DomainExt[1]) and
      (0.0 < coords[2] and coords[2] + size < DomainExt[2]);

  const bool isInCube =
      (coords[0] + size <= 0.0 or DomainExt[0] <= coords[0]) or
      (coords[1] + size <= 0.0 or DomainExt[1] <= coords[1]) or
      (coords[2] + size <= 0.0 or DomainExt[2] <= coords[2]);

  bool isOut = true;
  bool isIn = false;

  isOut &= isOutCube;
  isIn |= isInCube;

  constexpr int numSpheres = 3;
  for (int nSphere = 0; nSphere < numSpheres; nSphere++)
  {
    double clampCenter[3] = {centerSphere[nSphere][0],
                             centerSphere[nSphere][1],
                             centerSphere[nSphere][2]};

    for (int axis = 0; axis < DIM; ++axis)
      clampCenter[axis] = fmin(fmax(clampCenter[axis], coords[axis]), coords[axis] + size);

    double clampDist2 = 0.0;
    for (int axis = 0; axis < DIM; ++axis)
    {
      const double axisDist = clampCenter[axis] - centerSphere[nSphere][axis];
      clampDist2 += axisDist * axisDist;
    }

    const bool isOutSphere = (clampDist2 > RADIUS * RADIUS);

    bool isInSphere = true;
    for (int corner = 0; corner < (1u << DIM); ++corner)
    {
      double cornerDist2 = 0.0;
      for (int axis = 0; axis < DIM; ++axis)
      {
        const double axisDist = coords[axis] + bool(corner & (1u << axis)) * size - centerSphere[nSphere][axis];
        cornerDist2 += axisDist * axisDist;
      }

      isInSphere &= cornerDist2 <= RADIUS * RADIUS;
    }

    isOut &= isOutSphere;
    isIn |= isInSphere;
  }

  return (isOut ? ibm::OUT : isIn ? ibm::IN : ibm::INTERCEPTED);
}



/*
ibm::Partition DomainDecider(const double *physCoords, double physSize) {
  static bool XMortonOrder[nchild]{false, true, false, true, false, true, false, true};
  static bool YMortonOrder[nchild]{false, false, true, true, false, false, true, true};
  static bool ZMortonOrder[nchild]{false, false, false, false, true, true, true, true};
  double coords[3];
  std::array<bool, nchild> isOutsideDomain; // isOutside = true means that we want to retain the points
  isOutsideDomain.fill(true);
  // 1. At the initial stage the domain can be so skewed that it might lie within one element.
  // For that, those elements within which domain lies it returns as INTERCEPTED.
  // This solves the issues with strict and partial domain decider and we do not get 0 elements.
  if (creationStage == CREATION_STAGE::INITIAL) {
    coords[0] = physCoords[0] * scaling;
    coords[1] = physCoords[1] * scaling;
    coords[2] = physCoords[2] * scaling;
    if (((coords[0] <= 0.0) and (coords[0] + physSize * scaling >= DomainExt[0]))
        or ((coords[1] <= 0.0) and (coords[1] + physSize * scaling >= DomainExt[1]))
        or ((coords[2] <= 0.0) and (coords[2] + physSize * scaling >= DomainExt[2]))
        ) {
      return ibm::Partition::INTERCEPTED;
    }
  }
  // Check node by node basis
  for (int n = 0; n < nchild; n++) {
    coords[0] = (physCoords[0] + XMortonOrder[n] * physSize) * scaling;
    coords[1] = (physCoords[1] + YMortonOrder[n] * physSize) * scaling;
    coords[2] = (physCoords[2] + ZMortonOrder[n] * physSize) * scaling;
    // Check if inside any one of the sphere
    for (int nSphere = 0; nSphere < 3; nSphere++) {
      double dist = (centerSphere[nSphere][0] - coords[0]) * (centerSphere[nSphere][0] - coords[0]) +
                    (centerSphere[nSphere][1] - coords[1]) * (centerSphere[nSphere][1] - coords[1]) +
                    (centerSphere[nSphere][2] - coords[2]) * (centerSphere[nSphere][2] - coords[2]);
      if (dist < RADIUS * RADIUS) {
        isOutsideDomain[n] = false;
      }
    }
    // Check if inside the domain.
    // The domain extents are marked as IN boundary nodes.
    // It serves two purposes: 1. The element will automatically get marked Intercepted.
    // 2. The boundary nodes are IN nodes. So naturally transition to that.
    // Note that this the nomenclature is  opposite to that of spheres.
    // we want to retain the points that are inside domainBoundaries whereas in Sphere we
    // want to discard point that are inside the sphere.
    bool insideDomainBoundaries(true);
    for (int d = 0; d < DIM; d++) {
      insideDomainBoundaries = insideDomainBoundaries and (coords[d] > 0.0) and (coords[d] < DomainExt[d]);
    }
    isOutsideDomain[n] = isOutsideDomain[n] and (insideDomainBoundaries);
  }
  unsigned int numOutsidePoints = std::accumulate(isOutsideDomain.begin(), isOutsideDomain.end(), 0);
  if (numOutsidePoints == 0) { // No points is outside
    return ibm::IN;
  } else if (numOutsidePoints == nchild) { // All points are outside
    return ibm::OUT;
  }
  return ibm::INTERCEPTED;
}
*/

/**
 * main()
 */
int main(int argc, char *argv[]) {
  typedef unsigned int DENDRITE_UINT;
  PetscInitialize(&argc, &argv, NULL, NULL);
  DendroScopeBegin();
  _InitializeHcurve(DIM);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  dbgRank = rank;
  m_uiMaxDepth = 25;
  const DENDRITE_UINT eleOrder = 1;
  if (argc < 3) {
    if (not(rank)) {
      std::cout << "Usage: level maxLevel\n";
    }
    exit(EXIT_FAILURE);
  }
  const DENDRITE_UINT level = static_cast<DENDRITE_UINT>(std::atoi(argv[1]));
  maxLevel = static_cast<DENDRITE_UINT>(std::atoi(argv[2]));
  MPI_Comm comm = MPI_COMM_WORLD;
  using DTree = ot::DistTree<unsigned int, DIM>;
  DTree distTree = DTree::constructSubdomainDistTree(level, DomainDecider,
                                                     comm);
  /// ot::DA<DIM> *octDA = new ot::DA<DIM>(distTree, comm, eleOrder);
  std::cout << distTree.getTreePartFiltered().size() << "\n";

  if (rank == 0)
    printf("%*s[%d] Round %d, size %lu\n", rank, "", rank, 0, distTree.getTreePartFiltered().size());

  /// std::cout << octDA->getLocalElementSz() << "\n";
  creationStage = CREATION_STAGE::FINAL; /// Mark it as Final. Now we can not get 0 elements.
  for(int k = 0; k < 8; k++){
    std::vector<ot::OCT_FLAGS::Refine> refineFlags;

    /// generateRefinementFlags(octDA,distTree.getTreePartFiltered(),refineFlags);
    int local_refineAny = generateRefinementFlags(distTree.getTreePartFiltered(), refineFlags);

    int global_refineAny = 0;
    MPI_Allreduce(&local_refineAny, &global_refineAny, 1, MPI_INT, MPI_LOR, comm);

    if (global_refineAny)
    {
      ot::DistTree<unsigned int, DIM> newDistTree, surrDistTree;
      ot::DistTree<unsigned int, DIM>::distRemeshSubdomain(distTree, refineFlags, newDistTree, surrDistTree, ot::RemeshPartition::SurrogateInByOut, 0.3);
      /// ot::DA<DIM> *newDA = new ot::DA<DIM>(newDistTree, comm, eleOrder, 100, 0.3); //DistTree overload
      std::swap(distTree,newDistTree);
      /// std::swap(octDA,newDA);
      /// delete newDA;
      if(not(rank)) {
        /// std::cout << octDA->getLocalElementSz() << "\n";
        /// std::cout << distTree.getTreePartFiltered().size() << "\n";
      }
    }

    if (rank == 0)
      printf("%*s[%d] Round %d, size %lu\n", rank, "", rank, k+1, distTree.getTreePartFiltered().size());
  }

  ot::DA<DIM> *octDA = new ot::DA<DIM>(distTree, comm, eleOrder);

  if (rank == 0)
    printf("checkHangingNodes()\n");
  checkHangingNodes(octDA,distTree.getTreePartFiltered());

  /// // Huge ascii vtu files.
  /// io::vtk::oct2vtu(&(*distTree.getTreePartFiltered().cbegin()),
  ///                  distTree.getTreePartFiltered().size(),
  ///                  "octreeVtu",
  ///                  comm);

  DendroScopeEnd();
  PetscFinalize();
}
