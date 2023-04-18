#ifdef BUILD_WITH_PETSC
#include "petsc.h"
#endif
#include <iostream>
#include <distTree.h>
#include <oda.h>
#include <point.h>
#include <sfcTreeLoop_matvec_io.h>
#include <octUtils.h>
#include <filterFunction.h>
constexpr unsigned int DIM = 2;
constexpr unsigned int nchild = 1u << DIM;
static double xDomainExtent;
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
void printMaxCoords(ot::DA<DIM> & octDA, const std::vector<ot::TreeNode<unsigned, DIM>> &treePart)
{
  const size_t sz = octDA.getTotalNodalSz();
  auto partFront = octDA.getTreePartFront();
  auto partBack = octDA.getTreePartBack();
  const auto tnCoords = octDA.getTNCoords();
  {
    std::vector<double> maxCoords(4, 0.0);
    const int eleOrder = 1;
    const bool visitEmpty = false;
    const unsigned int padLevel = 0;
    const unsigned int npe = octDA.getNumNodesPerElement();
    ot::MatvecBaseCoords<DIM> loop(sz, eleOrder, visitEmpty, padLevel, tnCoords, &(*treePart.cbegin()), treePart.size(), *partFront, *partBack);
    while (!loop.isFinished())
    {
      if (loop.isPre() && loop.subtreeInfo().isLeaf())
      {
        const double *nodeCoordsFlat = loop.subtreeInfo().getNodeCoords();
        for (int i = 0; i < npe; i++)
        {
          for (int d = 0; d < DIM; d++)
          {
            if (maxCoords[d] < nodeCoordsFlat[i * DIM + d])
            {
              maxCoords[d] = nodeCoordsFlat[i * DIM + d];
            }
          }
        }
        loop.next();
      }
      else
      {
        loop.step();
      }
    }
    std::cout << "maxCoords == " << maxCoords[0] << " " << maxCoords[1] << " " << maxCoords[2] << "\n";
  }
}
void getBoundaryElements(const ot::DA<DIM>* octDA, const std::vector<ot::TreeNode<unsigned, DIM>> &treePart, unsigned  int eleOrder, const std::string fname){
  const size_t sz = octDA->getTotalNodalSz();
  auto partFront = octDA->getTreePartFront();
  auto partBack = octDA->getTreePartBack();
  const auto tnCoords = octDA->getTNCoords();
  int rank = octDA->getRankActive();
  int proc = octDA->getNpesActive();
  const unsigned int npe = octDA->getNumNodesPerElement();
  int counter = 0;
  std::string boundaryFname = fname+"_boundary" + std::to_string(rank) + "_" +std::to_string(proc) + ".txt";
  std::string nonBoundaryFname = fname+"_nonboundary"+std::to_string(rank) +"_" + std::to_string(proc) + ".txt";
  std::ofstream foutBoundary(boundaryFname.c_str());
  std::ofstream foutNonBoundary(nonBoundaryFname.c_str());
  ot::MatvecBaseCoords <DIM> loop(sz,eleOrder, false,0,tnCoords, &(*treePart.cbegin()), treePart.size(), *partFront,*partBack);
  while(!loop.isFinished()){
    if (loop.isPre() && loop.subtreeInfo().isLeaf()) {
      const double *nodeCoordsFlat = loop.subtreeInfo().getNodeCoords();
      if(loop.subtreeInfo().isElementBoundary()){
        counter++;
        for(int i = 0; i < npe; i++){
          foutBoundary << nodeCoordsFlat[2*i+0] << " "<< nodeCoordsFlat[2*i+1] <<  "\n";
        }
      }
      else{
        for(int i = 0; i < npe; i++) {
          foutNonBoundary << nodeCoordsFlat[2 * i + 0] << " " << nodeCoordsFlat[2 * i + 1] << "\n";
        }
      }
      loop.next();
    }
    else{
      loop.step();
    }
  }
  foutBoundary.close();
  foutNonBoundary.close();
  for(int i = 0; i < proc; i++){
    if(i == rank){
      std::cout << "Number of boundary elements in proc = " << i  << " = " << counter << "\n";
    }
    MPI_Barrier(octDA->getCommActive());
  }
}

ibm::Partition  DomainDecider(const double * p, double sz)
{
  const double boundsMin[2] = {0.0, 0.0};
  const double boundsMax[2] = {xDomainExtent, 1.0};

  const double cutoutMin[2] = {20.0/300.0, 20.0/300.0};
  const double cutoutMax[2] = {30.0/300.0, 30.0/300.0};

  using ibm::IN;
  using ibm::OUT;
  using ibm::INTERCEPTED;

  //
  // Keep out, discard in.
  //

  bool isOut =
  (
    ( p[0] > boundsMin[0] && p[0]+sz < boundsMax[0]  &&  p[1] > boundsMin[1] && p[1]+sz < boundsMax[1] )   // Inside bounds
      and
    ( p[0]+sz < cutoutMin[0] || p[0] > cutoutMax[0]  ||  p[1]+sz < cutoutMin[1] || p[1] > cutoutMax[1] )   // Outside cutout
  );


  bool isIn =
  (
    ( p[0]+sz <= boundsMin[0] || p[0] >= boundsMax[0]  ||  p[1]+sz <= boundsMin[1] || p[1] >= boundsMax[1] )  // Outside bounds
      or
    ( p[0] >= cutoutMin[0] && p[0]+sz <= cutoutMax[0]  &&  p[1] >= cutoutMin[1] && p[1]+sz <= cutoutMax[1] )  // Inside cutout
  );

  if (isIn && !isOut)
    return ibm::IN;
  else if (isOut && !isIn)
    return ibm::OUT;
  else if (!isIn && !isOut)
    return ibm::INTERCEPTED;
  else
    throw std::logic_error("Filter function returned both isIn and isOut.");
}


bool intersectsInnerBox(const double * p, double sz)
{
  const double cutoutMin[2] = {20.0/300.0, 20.0/300.0};
  const double cutoutMax[2] = {30.0/300.0, 30.0/300.0};

  return ( (p[0]+sz > cutoutMin[0] && p[0] < cutoutMax[0]) && (p[1]+sz > cutoutMin[1] && p[1] < cutoutMax[1]) );
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
  xDomainExtent = 0.5;
  const DENDRITE_UINT eleOrder = 1;
  if(argc < 3){
    if(not(rank)){
      std::cout << "Usage: level xExtent\n";
    }
    exit(EXIT_FAILURE);
  }
  const DENDRITE_UINT level = static_cast<DENDRITE_UINT>(std::atoi(argv[1]));
  xDomainExtent =  static_cast<double>(std::atof(argv[2]));
  std::cout << level << " "<<  xDomainExtent << "\n";
  MPI_Comm comm = MPI_COMM_WORLD;
  constexpr bool printTreeOn = false;  // Can print the contents of the tree vectors.
  unsigned int extents[] = {1,2,1};
  std::array<unsigned int,DIM> a;
  for (int d = 0; d < DIM; ++d)
    a[d] = extents[d];
  using DTree = ot::DistTree<unsigned int, DIM>;
  DTree distTree = DTree::constructSubdomainDistTree( level,DomainDecider,
                                                      comm);
  ot::DA<DIM> *octDA = new ot::DA<DIM>(distTree, comm, eleOrder);
  /// printMaxCoords(*octDA, distTree.getTreePartFiltered());
  size_t oldTreeSize = 0;
  size_t refinedTreeSize = 0;
  // Access the original tree as a list of tree nodes.
  {
    const std::vector<ot::TreeNode<unsigned int, DIM>> &treePart = distTree.getTreePartFiltered();
    oldTreeSize = treePart.size();
    /// if (printTreeOn)
    ///   printTree(treePart, level+1);
  }
  std::cout << "Old Tree \n";
  std::cout << "Num elements: " << oldTreeSize << "\n";
  getBoundaryElements(octDA, distTree.getTreePartFiltered(),eleOrder, "Old");
  std::string fname = "BoNodes"+std::to_string(rank) + ".txt";
  std::ofstream fout(fname.c_str());
//  DomainInfo fullDomain;
//  fullDomain.min.fill(0.0);
//  fullDomain.max.fill(1.0);
//  IO::writeBoundaryElements(octDA,"subDA",fullDomain);
for(int i = 0; i < 7; i++) {
  std::vector<ot::OCT_FLAGS::Refine> refineFlags(octDA->getLocalElementSz(), ot::OCT_FLAGS::Refine::OCT_NO_CHANGE);
  const size_t sz = octDA->getTotalNodalSz();
  auto partFront = octDA->getTreePartFront();
  auto partBack = octDA->getTreePartBack();
  const auto tnCoords = octDA->getTNCoords();
  ot::MatvecBaseCoords<DIM> loop(sz, octDA->getElementOrder(), false, 0, tnCoords, &(*distTree.getTreePartFiltered().cbegin()), distTree.getTreePartFiltered().size(), *partFront, *partBack);
  int counter = 0;
  while (!loop.isFinished()) {
    if (loop.isPre() && loop.subtreeInfo().isLeaf()) {

      double physCoords[DIM];
      double physSize;
      treeNode2Physical(loop.subtreeInfo().getCurrentSubtree(), physCoords, physSize);
      if (intersectsInnerBox(physCoords, physSize))
        refineFlags[counter] = ot::OCT_FLAGS::Refine::OCT_REFINE;

      /// const double *nodeCoordsFlat = loop.subtreeInfo().getNodeCoords();
      /// for (int numEle = 0; numEle < 4; numEle++) {
      ///   if ((nodeCoordsFlat[numEle * DIM + 0] >= (20. / 300.)) and
      ///       (nodeCoordsFlat[numEle * DIM + 0] <= (30. / 300.)) and
      ///       (nodeCoordsFlat[numEle * DIM + 1] >= (20. / 300.)) and
      ///       (nodeCoordsFlat[numEle * DIM + 1] <= (30. / 300.))) {
      ///     if (loop.subtreeInfo().isElementBoundary()) {
      ///       refineFlags[counter] = ot::OCT_FLAGS::Refine::OCT_REFINE;
      ///     }
      ///   }
      /// }

      counter++;
      loop.next();
    } else {
      loop.step();
    }
  }
  
  if (!(counter == octDA->getLocalElementSz()))
  {
    printf("WARNING!!  counter==%d, octDA->getLocalElementSz()==%d\n", counter, (int) octDA->getLocalElementSz());
  }


  // distRemeshSubdomain()
  ot::DistTree<unsigned int, DIM> newDistTree, surrDistTree;
  ot::DistTree<unsigned int, DIM>::distRemeshSubdomain(distTree, refineFlags, newDistTree, surrDistTree, ot::RemeshPartition::SurrogateInByOut, 0.3);
  // =======================================================================
  // Currently the boundary is detected from applying the carving function.
  // The boundary is not detected based solely on the existence of TreeNodes
  // in the tree. So DA() constructor needs DistTree, i.e.,
  // newDistTree instead of newTree. Otherwise, DA() constructor uses
  // default boundary definition of the unit cube.
  // =======================================================================
  ot::DA<DIM> *newDA = new ot::DA<DIM>(newDistTree, comm, eleOrder, 100, 0.3); //DistTree overload
  /// printMaxCoords(*newDA, newDistTree.getTreePartFiltered());
  // Access the refined tree as a list of tree nodes.
  {
    const std::vector<ot::TreeNode<DENDRITE_UINT, DIM>> &newTree = newDistTree.getTreePartFiltered();
    /// const std::vector<ot::TreeNode<DENDRITE_UINT, DIM>> &surrTree = surrDistTree.getTreePartFiltered();
    refinedTreeSize = newTree.size();
    /// if (printTreeOn)
    ///   printTree(newTree, level+1);
  }
  std::swap(octDA, newDA);
  delete newDA;
  std::swap(distTree,newDistTree);
//  IO::writeBoundaryElements(octDA, ("RefinedsubDA"+std::to_string(i)).c_str(), fullDomain);
  std::cout << "New tree \n";
  std::cout << "Num elements: " << refinedTreeSize << "\n";
  getBoundaryElements(octDA, distTree.getTreePartFiltered(), eleOrder, ("Refined"+std::to_string(i)).c_str());
}
  delete octDA;
  _DestroyHcurve();
  DendroScopeEnd();
  PetscFinalize();
}
