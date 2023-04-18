//
// Created by maksbh on 7/2/20.
//

#include <oda.h>
#define DIM 3
int main(int argc, char * argv[]){

  PetscInitialize(&argc, &argv, NULL, NULL);
  _InitializeHcurve(DIM);
  /** initial setup of DA **/
  int rank; MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  int eleOrder = 1;
  m_uiMaxDepth = 10;
  int level = 4;
  std::vector<ot::TreeNode<unsigned int, DIM>> treePart;
  ot::createRegularOctree(treePart, level, MPI_COMM_WORLD);
  ot::DA<DIM> *octDA = new ot::DA<DIM>(treePart, MPI_COMM_WORLD, eleOrder);


  /** Point whose processor needs to be determined **/

  const double point[3]{0.2,0.2,0.2};

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
  treeNode[0] = ot::TreeNode<unsigned int,DIM>(treeNodeArray,m_uiMaxDepth);

  /** Always getting 0 as the owner proc**/
  int ownerRank[1];
  octDA->computeTreeNodeOwnerProc(treeNode.data(),1,ownerRank);

  std::cout << ownerRank[0] << "\n";

  delete octDA;
  PetscFinalize();
}