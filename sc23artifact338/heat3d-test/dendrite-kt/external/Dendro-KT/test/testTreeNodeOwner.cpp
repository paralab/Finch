

#ifdef BUILD_WITH_PETSC
#include "petsc.h"
#endif

#include <oda.h>

#define DIM 3

int main(int argc, char * argv[])
{
  PetscInitialize(&argc, &argv, NULL, NULL);
  DendroScopeBegin();
  _InitializeHcurve(DIM);

  MPI_Comm comm = MPI_COMM_WORLD;

  /** initial setup of DA **/
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  int eleOrder = 1;
  m_uiMaxDepth = 10;
  int level = 4;

  std::vector<ot::TreeNode<unsigned int, DIM>> treePart;
  ot::createRegularOctree(treePart, level, MPI_COMM_WORLD);

  /** Print out the number of elements on each processor. **/
  MPI_Barrier(comm);
  for (int r = 0; r < size; ++r)
  {
    if (r == rank)
      printf("Rank [%d]: #elements==%lu\n", rank, (unsigned long) treePart.size());
    MPI_Barrier(comm);
  }

  ot::DA<DIM> *octDA = new ot::DA<DIM>(ot::DistTree<unsigned, DIM>(treePart, MPI_COMM_WORLD), MPI_COMM_WORLD, eleOrder);

  /** Print out the number of nodes on each processor. **/
  MPI_Barrier(comm);
  for (int r = 0; r < size; ++r)
  {
    if (r == rank)
      printf("Rank [%d]: #nodes==%lu\n", rank, (unsigned long) octDA->getLocalNodalSz());
    MPI_Barrier(comm);
  }


  /** Point whose processor needs to be determined **/
  const double point[3]{0.95,0.95,0.95}; 

  /** Assume the domain is (0,0,0) to  (1,1,1). Converting from physical to octant coordinates **/
  const double physToOct[DIM] = {
      (1u << (m_uiMaxDepth)) / 1.0,
      (1u << (m_uiMaxDepth)) / 1.0,
      (1u << (m_uiMaxDepth)) / 1.0,
  };
  unsigned int octCoords[DIM];
  for (unsigned int dim = 0; dim < DIM; dim++) {
    octCoords[dim] = (unsigned int) (point[dim] * physToOct[dim]);
  }

  /*** Creating tree Node **/
  std::array<unsigned int,DIM> treeNodeArray;
  for(unsigned int dim = 0; dim <DIM; dim++){
    treeNodeArray[dim] = octCoords[dim];
  }

  std::vector<ot::TreeNode<unsigned int,DIM>> treeNode(1);

  /** The logic is that we want to compare only on based on coordinates. So, we put the tree Node to maxDepth**/
  treeNode[0] = ot::TreeNode<unsigned int,DIM>(treeNodeArray, m_uiMaxDepth);

  /** Always getting 0 as the owner proc**/
  int ownerRank[1];
  octDA->computeTreeNodeOwnerProc(treeNode.data(), 1, ownerRank);
  printf("[%d]: ownerRank[0]==%d\n", rank, ownerRank[0]);

  delete octDA;
  _DestroyHcurve();
  DendroScopeEnd();
  PetscFinalize();
}
