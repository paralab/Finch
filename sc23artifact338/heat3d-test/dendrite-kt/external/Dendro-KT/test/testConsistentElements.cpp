#ifdef BUILD_WITH_PETSC
#include "petsc.h"
#endif
#include <iostream>
#include <distTree.h>
#include <oda.h>
#include <point.h>
#include <sfcTreeLoop_matvec_io.h>
#include <octUtils.h>
constexpr unsigned int DIM = 2;
constexpr unsigned int nchild = 1u << DIM;
template <unsigned int dim>
void printTree(const std::vector<ot::TreeNode<unsigned int, dim>> &treePart, int elev)
{
  std::cout << "Tree\n";
  for (const ot::TreeNode<unsigned int, dim> &tn : treePart)
  {
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
void generateRefinementFlags(const ot::DA<DIM> * octDA, const std::vector<ot::TreeNode<unsigned, DIM>> &treePart, std::vector<ot::OCT_FLAGS::Refine> & refineFlags){
  const size_t sz = octDA->getTotalNodalSz();
  auto partFront = octDA->getTreePartFront();
  auto partBack = octDA->getTreePartBack();
  const auto tnCoords = octDA->getTNCoords();
  const unsigned int eleOrder = octDA->getElementOrder();
  refineFlags.resize(octDA->getLocalElementSz());
  const unsigned int npe = octDA->getNumNodesPerElement();
  int counter = 0;
  ot::MatvecBaseCoords <DIM> loop(sz,eleOrder, false,0,tnCoords,  &(*treePart.cbegin()), treePart.size(), *partFront,*partBack);
  while(!loop.isFinished()){
    if (loop.isPre() && loop.subtreeInfo().isLeaf()) {
      if (loop.subtreeInfo().isElementBoundary()) {
        refineFlags[counter] = ot::OCT_FLAGS::Refine::OCT_REFINE;
      } else {
        refineFlags[counter] = ot::OCT_FLAGS::Refine::OCT_NO_CHANGE;
      }
      counter++;
      loop.next();
    }
    else{
      loop.step();
    }
  }
}
unsigned int getNumElements(const ot::DA<DIM> * octDA, const std::vector<ot::TreeNode<unsigned, DIM>> &treePart){
  const size_t sz = octDA->getTotalNodalSz();
  auto partFront = octDA->getTreePartFront();
  auto partBack = octDA->getTreePartBack();
  const auto tnCoords = octDA->getTNCoords();
  const unsigned int eleOrder = octDA->getElementOrder();
  const unsigned int npe = octDA->getNumNodesPerElement();
  int counter = 0;
  ot::MatvecBaseCoords <DIM> loop(sz,eleOrder, false,0,tnCoords, &(*treePart.cbegin()), treePart.size(),  *partFront,*partBack);
  while(!loop.isFinished()){
    if (loop.isPre() && loop.subtreeInfo().isLeaf()) {
      counter++;
      loop.next();
    }
    else{
      loop.step();
    }
  }
  return counter;
}
ibm::Partition DomainDecider(const double * physCoords, double physSize)
{
  static bool XMortonOrder[4] {false, true ,false, true};
  static bool YMortonOrder[4] {false, false ,true, true};
  double coords[2];
  static const unsigned int numChildren = 4;
  std::array<bool, numChildren> isOutsideDomain; // isOutside = true means that we want to retain the points
  double centerCircle[2] {0.0625,0.125};
  double radiusCircle = 0.0625;
  isOutsideDomain.fill(true);
  for(int n = 0; n < numChildren; n++) {
    coords[0] = physCoords[0] + XMortonOrder[n] * physSize;
    coords[1] = physCoords[1] + YMortonOrder[n] * physSize;
    // Check if inside the domain. The domain is [0,1] X [0,0.25]
    // The domain extents are marked as IN boundary nodes.
    // It serves two purposes: 1. The element will automatically get marked Intercepted.
    // 2. The boundary nodes are IN nodes. So naturally transition to that.
    // This way has issues with strict and partial domain decider but we can worry about it later.
    if ((coords[1] >= 0.25) or (coords[0] == 0) or (coords[1] == 0) or (coords[0] == 1 )) {
      isOutsideDomain[n] = false;
    }
    // Check if inside the circle
    double dist = (centerCircle[0] - coords[0]) * (centerCircle[0] - coords[0]) +
                  (centerCircle[1] - coords[1]) * (centerCircle[1] - coords[1]);
    if (dist < radiusCircle*radiusCircle) {
      isOutsideDomain[n] = false;
    }
  }
  unsigned int numOutsidePoints = std::accumulate(isOutsideDomain.begin(),isOutsideDomain.end(),0);
  if(numOutsidePoints == 0){ // No points is outside
    return ibm::IN;
  }
  else if(numOutsidePoints == 4){ // All points are outside
    return ibm::OUT;
  }
  return ibm::INTERCEPTED;
}
/**
 * main()
 */
int main(int argc, char * argv[]){
  typedef unsigned int DENDRITE_UINT;
  PetscInitialize(&argc, &argv, NULL, NULL);
  DendroScopeBegin();
  _InitializeHcurve(DIM);
  int rank; MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  m_uiMaxDepth = 10;
  const DENDRITE_UINT eleOrder = 1;
  if(argc < 2){
    if(not(rank)){
      std::cout << "Usage: level\n";
    }
    exit(EXIT_FAILURE);
  }
  const DENDRITE_UINT level = static_cast<DENDRITE_UINT>(std::atoi(argv[1]));
  MPI_Comm comm = MPI_COMM_WORLD;
  using DTree = ot::DistTree<unsigned int, DIM>;
  DTree distTree = DTree::constructSubdomainDistTree( level,DomainDecider,
                                                      comm);
  const std::vector<ot::TreeNode<unsigned int,DIM>> & treePart = distTree.getTreePartFiltered();

  ot::DA<DIM> *octDA = new ot::DA<DIM>(distTree, comm, eleOrder);  //TODO in release mode this segfaults, why?

  std::cout << "da elements: " << octDA->getLocalElementSz()
            << "  loop elements: "<< getNumElements(octDA, treePart)  << "  "
            << (octDA->getLocalElementSz() == getNumElements(octDA, treePart) ? "(equal)" : "(not equal!)") << "\n";
  assert(octDA->getLocalElementSz() == getNumElements(octDA, treePart));

  std::vector<ot::OCT_FLAGS::Refine> refineFlags;
  generateRefinementFlags(octDA, treePart,refineFlags);
  ot::DistTree<unsigned int, DIM> newDistTree, surrDistTree;
  ot::DistTree<unsigned int, DIM>::distRemeshSubdomain(distTree, refineFlags, newDistTree, surrDistTree, ot::RemeshPartition::SurrogateInByOut, 0.3);
  const std::vector<ot::TreeNode<unsigned int, DIM>> &newTreePart = newDistTree.getTreePartFiltered();
  ot::DA<DIM> *newDA = new ot::DA<DIM>(newDistTree, comm, eleOrder, 100, 0.3); //DistTree overload

  std::cout << "da elements: " << newDA->getLocalElementSz()
            << "  loop elements: "<< getNumElements(newDA, newTreePart) <<  "  "
            << (newDA->getLocalElementSz() == getNumElements(newDA, newTreePart) ? "(equal)" : "(not equal!)") << "\n";
  assert(newDA->getLocalElementSz() == getNumElements(newDA, newTreePart)); // assertion fails using original sfctreeloop
  DendroScopeEnd();
  PetscFinalize();
}
