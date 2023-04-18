#include <iostream>
#include <distTree.h>
#include <oda.h>
#include <point.h>
#include <sfcTreeLoop_matvec_io.h>

#include <octUtils.h>


constexpr unsigned int DIM = 2;
constexpr unsigned int nchild = 1 << DIM;

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



void printMaxCoords(ot::DA<DIM> & octDA, const std::vector<ot::TreeNode<unsigned, DIM>> &treePart){
  const size_t sz = octDA.getTotalNodalSz();
  auto partFront = octDA.getTreePartFront();
  auto partBack = octDA.getTreePartBack();
  const auto tnCoords = octDA.getTNCoords();
  {
    std::vector<double> maxCoords(4, 0.0);
    ot::MatvecBaseCoords<DIM> loop(sz, 1, false, 0, tnCoords, &(*treePart.cbegin()), treePart.size(), *partFront, *partBack);
    while (!loop.isFinished()) {
      if (loop.isPre() && loop.subtreeInfo().isLeaf()) {
        const double *nodeCoordsFlat = loop.subtreeInfo().getNodeCoords();
        for (int i = 0; i < nchild; i++) {
          for (int d = 0; d < DIM; d++) {
            if (maxCoords[d] < nodeCoordsFlat[i * DIM + d]) {
              maxCoords[d] = nodeCoordsFlat[i * DIM + d];
            }
          }
        }
        loop.next();
      } else {
        loop.step();
      }
    }
    std::cout << "maxCoords == " << maxCoords[0] << " " << maxCoords[1] << " " << maxCoords[2] << "\n";
  }
}



int main(int argc, char * argv[]){
  PetscInitialize(&argc, &argv, NULL, NULL);
  DendroScopeBegin();
  typedef unsigned int DENDRITE_UINT;
  _InitializeHcurve(DIM);
  int eleOrder = 1;
  int ndof = 1;
  m_uiMaxDepth = 10;
  int level = 3;
  std::vector<ot::TreeNode<unsigned int, DIM>> treePart;

  constexpr bool printTreeOn = false;  // Can print the contents of the tree vectors.

  unsigned int extents[] = {1,2,1};
  std::array<double,DIM> a;
  for (int d = 0; d < DIM; ++d)
    a[d] = double(1u << extents[d]) / double(1u << level);

  using DTree = ot::DistTree<unsigned int, DIM>;
  DTree distTree = DTree::constructSubdomainDistTree( level,
                                                      DTree::BoxDecider(a),
                                                      MPI_COMM_WORLD);

  treePart = distTree.getTreePartFiltered();
  ot::DA<DIM> *octDA = new ot::DA<DIM>(distTree, MPI_COMM_WORLD, eleOrder);

  //Old way, doesn't support refining subda.
  /// ot::DA<DIM> *octDA = new ot::DA<DIM>();
  /// ot::constructRegularSubdomainDA<DIM>(*octDA,treePart,level,a,eleOrder,MPI_COMM_WORLD);

  if (printTreeOn)
    printTree(treePart, level+1);
  printMaxCoords(*octDA, distTree.getTreePartFiltered());

  std::vector<ot::OCT_FLAGS::Refine> refineFlags(treePart.size(),ot::OCT_FLAGS::Refine::OCT_REFINE);

  // distRemeshSubdomain()
  ot::DistTree<unsigned int, DIM> newDistTree, surrDistTree;
  ot::DistTree<unsigned int, DIM>::distRemeshSubdomain(distTree, refineFlags, newDistTree, surrDistTree, ot::RemeshPartition::SurrogateInByOut, 0.3);
  const std::vector<ot::TreeNode<DENDRITE_UINT, DIM>> &newTree = newDistTree.getTreePartFiltered();
  const std::vector<ot::TreeNode<DENDRITE_UINT, DIM>> &surrTree = surrDistTree.getTreePartFiltered();

  //Old way, returns back the whole domain.
  /// std::vector<ot::TreeNode<DENDRITE_UINT, DIM>>  newTree;
  /// std::vector<ot::TreeNode<DENDRITE_UINT, DIM>>  surrTree;
  /// ot::SFC_Tree<DENDRITE_UINT , DIM>::distRemeshWholeDomain(treePart, refineFlags, newTree, surrTree, 0.3, ot::RemeshPartition::SurrogateInByOut, MPI_COMM_WORLD);

  std::cout << "\n-------\n";

  if (printTreeOn)
    printTree(newTree, level+1);
  ot::DA<DIM> * newDA =new ot::DA<DIM>(newDistTree,MPI_COMM_WORLD,eleOrder,100,0.3);
  printMaxCoords(*newDA, newDistTree.getTreePartFiltered());

  std::swap(octDA, newDA);
  delete newDA;
  DendroScopeEnd();
  PetscFinalize();
}
