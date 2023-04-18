//
// Created by maksbh on 9/16/19.
//

#include <iostream>
#include <oda.h>
#include <point.h>

#define DIM 2
int main(int argc, char * argv[]){

  PetscInitialize(&argc, &argv, NULL, NULL);
  _InitializeHcurve(DIM);
  int rank; MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  int eleOrder = 1;
  m_uiMaxDepth = 10;
  int level = 4;
  std::vector<ot::TreeNode<unsigned int, DIM>> treePart;
  ot::createRegularOctree(treePart, level, MPI_COMM_WORLD);
  ot::DA<DIM> *octDA = new ot::DA<DIM>(treePart, MPI_COMM_WORLD, eleOrder);

  std::cout << octDA->getLocalNodalSz() << "\n";

  PetscFinalize();
}
