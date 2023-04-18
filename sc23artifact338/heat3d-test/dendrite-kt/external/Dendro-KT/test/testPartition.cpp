// @author maksbh
// @note mpirun -np 16 $exec 6 8
// @note The code will output following file
//   output0.txt  output1.txt  outputInitial.txt  outputWeighted.txt
// with the outputInitial load distribution of initial mesh.
// output0 and output1 with load distribution after 1 and 2 iter of refinement
// and outputWeighted after the weighted partitioning
// @note Additionally, the code will output
//   BoundaryCoords_0.txt  BoundaryCoords_1.txt  BoundaryCoordsInital.txt  BoundaryCoordsWighted.txt
//  with the boundary coords.
// @note gnuplot command: splot 'BoundaryCoordsWighted.txt' u 1:2:3 w p 

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
//#include <IO/VTU.h>
constexpr unsigned int DIM = 3;
typedef ot::TreeNode<unsigned int, DIM> TREENODE;
constexpr unsigned int nchild = 1u << DIM;
static constexpr double RADIUS = 0.45;
static constexpr double scaling = 1.0; /// The octant Coords are scaled by 16 in each direction.
static constexpr std::array<double, 3> DOMAIN{1.0/16, 1.0/16, 4.0/16}; /// Physical domain of 1 X 1 X 16 is carved out from [1,1,16]
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
  refineFlags.resize(octDA->getLocalElementSz());
  std::fill(refineFlags.begin(), refineFlags.end(), ot::OCT_FLAGS::Refine::OCT_REFINE);
}
static void PrintBoundaryNodes(ot::DA<DIM> * octDA, const std::string & filename){
  std::vector<std::size_t> bdy_index;
  octDA->getBoundaryNodeIndices(bdy_index);
  double coords[DIM];
  if(octDA->getRankAll() == 0){
    std::ofstream fout(filename);
    fout.close();
  }
  MPI_Barrier(MPI_COMM_WORLD);
  for(int rank = 0; rank < octDA->getNpesAll(); rank++){
    if(rank == octDA->getRankAll()){
      std::ofstream fout(filename,std::ios::app);
      for (int i = 0; i < bdy_index.size(); i++) {
        ot::treeNode2Physical(octDA->getTNCoords()[bdy_index[i] + octDA->getLocalNodeBegin()],octDA->getElementOrder(),coords);
        fout << coords[0] << " " << coords[1] << " " << coords[2] << "\n";
      }
      fout.close();
      MPI_Barrier(MPI_COMM_WORLD);
    }
  }
}
static const auto DomainDecider = ot::DistTree<unsigned, 3>::BoxDecider(DOMAIN);

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
  PrintBoundaryNodes(octDA,"BoundaryCoordsInital.txt");
  creationStage = CREATION_STAGE::FINAL; /// Mark it as Final. Now we can not get 0 elements.
  for(int k = 0; k < 2; k++){
    std::vector<ot::OCT_FLAGS::Refine> refineFlags;
    generateRefinementFlags(octDA,treePart,refineFlags);
    ot::DistTree<unsigned int, DIM> newDistTree, surrDistTree;
    ot::DistTree<unsigned int, DIM>::distRemeshSubdomain(distTree, refineFlags, newDistTree, surrDistTree, ot::RemeshPartition::SurrogateInByOut, 0.3);
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
    MPI_Gather(&lsz,1,MPI_INT,lelem.data(),1,MPI_INT,0,comm);
    if(not(rank)) {
      std::ofstream fout("output"+ std::to_string(k) +".txt" );
      for (int i = 0; i < lelem.size(); i++) {
        fout << i << " " << lelem[i] << "\n";
      }
      fout.close();
    }
    PrintBoundaryNodes(octDA,"BoundaryCoords_" + std::to_string(k) + ".txt");
  }
  /// Weighted partitioning
  auto tree = distTree.getTreePartFiltered();
  par::partitionW(tree, par::defaultWeight,comm);
  DTree newDtree(tree, MPI_COMM_WORLD);
  newDtree.filterTree(distTree.getDomainDecider());
  delete octDA;
  ot::DA<DIM> *newDA = new ot::DA<DIM>(newDtree, comm, eleOrder, 100, 0.3); //DistTree overload
  std::swap(octDA,newDA);
  lsz = octDA->getLocalElementSz();
  MPI_Gather(&lsz,1,MPI_INT,lelem.data(),1,MPI_INT,0,comm);
  if(not(rank)) {
    std::ofstream fout("outputWeighted.txt" );
    for (int i = 0; i < lelem.size(); i++) {
      fout << i << " " << lelem[i] << "\n";
    }
    fout.close();
  }
  PrintBoundaryNodes(octDA,"BoundaryCoordsWighted.txt");
  DendroScopeEnd();
  PetscFinalize();
}
