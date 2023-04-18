//#ifdef BUILD_WITH_PETSC
//
//#include "petsc.h"
//
//#endif

#include <iostream>
#include <distTree.h>
#include <oda.h>
#include <point.h>
#include <sfcTreeLoop_matvec_io.h>
#include <octUtils.h>
#include <parUtils.h>
#include <IO/VTU.h>


//constexpr unsigned int DIM = 3;
typedef ot::TreeNode<unsigned int, DIM> TREENODE;
constexpr unsigned int nchild = 1u << DIM;
static constexpr double RADIUS = 0.45;

static constexpr double scaling = 16.0; /// The octant Coords are scaled by 16 in each direction.
/// The three spheres in the domain
static double centerSphere[3][3]{{0.5, 0.5, 8.0},
                                 {0.5, 0.0, 8.0 + RADIUS * 2 + 0.01},
                                 {0.5, 1.0, 8.0 + RADIUS * 2 + 0.01}};

static constexpr double DOMAIN[3]{1.0, 1.0, 4.0}; /// Physical domain of 1 X 1 X 16 is carved out from [1,1,16]
static unsigned int maxLevel;


enum CREATION_STAGE : bool {
  INITIAL = true,
  FINAL = false
};
static CREATION_STAGE creationStage = CREATION_STAGE::INITIAL;


/**
 * This function only refines the elements marked as boundary octants
 * @param octDA
 * @param refineFlags
 */
void generateRefinementFlags(const ot::DA<DIM> *octDA,const std::vector<TREENODE> &treePart,
                             std::vector<ot::OCT_FLAGS::Refine> &refineFlags) {
  const size_t sz = octDA->getTotalNodalSz();
  auto partFront = octDA->getTreePartFront();
  auto partBack = octDA->getTreePartBack();
  const auto tnCoords = octDA->getTNCoords();
  const unsigned int eleOrder = octDA->getElementOrder();
  refineFlags.resize(octDA->getLocalElementSz());
  std::fill(refineFlags.begin(), refineFlags.end(), ot::OCT_FLAGS::Refine::OCT_NO_CHANGE);
  const unsigned int npe = octDA->getNumNodesPerElement();
  int counter = 0;

  ot::MatvecBaseCoords<DIM> loop(sz, eleOrder, false, 0, tnCoords, &(*treePart.cbegin()), treePart.size(), *partFront,*partBack);
  while (!loop.isFinished()) {
    if (loop.isPre() && loop.subtreeInfo().isLeaf()) {
      refineFlags[counter] = ot::OCT_FLAGS::Refine::OCT_REFINE;
      counter++;
      loop.next();

    } else {
      loop.step();
    }
  }
}



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
    if (((coords[0] <= 0.0) and (coords[0] + physSize * scaling >= DOMAIN[0]))
        or ((coords[1] <= 0.0) and (coords[1] + physSize * scaling >= DOMAIN[1]))
        or ((coords[2] <= 0.0) and (coords[2] + physSize * scaling >= DOMAIN[2]))
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
      insideDomainBoundaries = insideDomainBoundaries and (coords[d] > 0.0) and (coords[d] < DOMAIN[d]);
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


/**
 * main()
 */
int main(int argc, char *argv[]) {
  typedef unsigned int DENDRITE_UINT;
  PetscInitialize(&argc, &argv, NULL, NULL);
  _InitializeHcurve(DIM);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
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
  const std::vector<ot::TreeNode<unsigned int, DIM>> &treePart = distTree.getTreePartFiltered();


  ot::DA<DIM> *octDA = new ot::DA<DIM>(distTree, comm, eleOrder);
  std::vector<int> lelem(octDA->getNpesAll());
  int lsz = octDA->getLocalElementSz();
  MPI_Gather(&lsz,1,MPI_INT,lelem.data(),1,MPI_INT,0,comm);
  if(not(rank)) {
    std::ofstream fout("outputInitial.txt");
    for (int i = 0; i < lelem.size(); i++) {
      fout << i << " " << lelem[i] << "\n";
    }
    fout.close();
  }
  MPI_Barrier(MPI_COMM_WORLD);

  creationStage = CREATION_STAGE::FINAL; /// Mark it as Final. Now we can not get 0 elements.

  for(int k = 0; k < 2; k++){
    std::vector<ot::OCT_FLAGS::Refine> refineFlags;
    generateRefinementFlags(octDA,treePart,refineFlags);
    ot::DistTree<unsigned int, DIM> newDistTree, surrDistTree;
    ot::DistTree<unsigned int, DIM>::distRemeshSubdomain(distTree, refineFlags, newDistTree, surrDistTree,ot::SurrogateInByOut, 0.3);


    ot::DA<DIM> *newDA = new ot::DA<DIM>(newDistTree, comm, eleOrder, 100, 0.3); //DistTree overload
    std::swap(distTree,newDistTree);
    std::swap(octDA,newDA);
    delete newDA;
    int lsz = octDA->getLocalElementSz();
    for (int i = 0; i < octDA->getNpesAll(); i++) {
      if(octDA->getRankAll() == i){
        std::cout << k << " " << i << " " << octDA->getLocalElementSz() << "\n";
      }
      MPI_Barrier(comm);
    }
    std::cout << "Finished k\n";
    MPI_Gather(&lsz,1,MPI_INT,lelem.data(),1,MPI_INT,0,comm);
    if(not(rank)) {
      std::ofstream fout("output"+ std::to_string(k) +".txt" );
      for (int i = 0; i < lelem.size(); i++) {

        fout << i << " " << lelem[i] << "\n";
      }
      fout.close();
    }

  }
  auto tree = distTree.getTreePartFiltered();
  par::partitionW(tree, par::defaultWeight,comm);
  DTree newDtree(tree, MPI_COMM_WORLD);

  ot::DA<DIM> *newDA = new ot::DA<DIM>(newDtree, comm, eleOrder, 100, 0.3); //DistTree overload
  std::swap(octDA,newDA);
  lsz = octDA->getLocalElementSz();
  MPI_Gather(&lsz,1,MPI_INT,lelem.data(),1,MPI_INT,0,comm);
  if(not(rank)) {
    std::ofstream fout("outputFinal.txt" );
    for (int i = 0; i < lelem.size(); i++) {

      fout << i << " " << lelem[i] << "\n";
    }
    fout.close();
  }
//  checkHangingNodes(octDA,distTree.getTreePartFiltered());
  std::array<DENDRITE_REAL,DIM> min; /// Minimum of the domain
  std::array<DENDRITE_REAL,DIM> max; /// Minimum of the domain
  min.fill(0.0);
  max.fill(16.0);
  DomainInfo fullDADomain{min,max};
  DomainExtents extents(fullDADomain);
  IO::writeBoundaryElements(octDA,newDtree.getTreePartFiltered(),"boundary","bnd",extents);
  PetscFinalize();

}


