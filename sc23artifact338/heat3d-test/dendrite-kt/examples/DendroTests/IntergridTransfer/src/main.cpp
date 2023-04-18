//
// Created by maksbh on 11/9/20.
//

#include <iostream>
#include <oda.h>
#include <point.h>
#include <sfcTreeLoop_matvec_io.h>
#include <intergridTransfer.h>
static constexpr unsigned int DIM = 2;
typedef ot::DA<DIM> DA;
typedef unsigned int DENDRITE_UINT;

typedef ot::DistTree<unsigned int, DIM> DistTREE;
typedef ot::TreeNode<DENDRITE_UINT, DIM> TREENODE;

/**
 * @brief Generate Flags that refine only the boundary elements.
 * If you set all to Refine. The code works.
 * @param octDA
 * @param refineFlags
 */
void generateRefinementFlags(ot::DA<DIM> * octDA, std::vector<ot::OCT_FLAGS::Refine> & refineFlags,const ot::DistTree<unsigned int, DIM> &distTree){
    const size_t sz = octDA->getTotalNodalSz();
    auto partFront = octDA->getTreePartFront();
    auto partBack = octDA->getTreePartBack();
    const auto tnCoords = octDA->getTNCoords();
    ot::MatvecBaseCoords <DIM> loop(sz,octDA->getElementOrder(), false,0,tnCoords,&(*distTree.getTreePartFiltered().cbegin()), distTree.getTreePartFiltered().size(), *partFront,*partBack);
    int counter = 0;
    while(!loop.isFinished()){
        if (loop.isPre() && loop.subtreeInfo().isLeaf()) {
            bool boundaryOctant = loop.subtreeInfo().isElementBoundary();
            if(boundaryOctant){
                refineFlags[counter]=ot::OCT_FLAGS::Refine::OCT_REFINE;
            }
            else{
                refineFlags[counter]=ot::OCT_FLAGS::Refine::OCT_NO_CHANGE;
            }
            counter++;
            loop.next();
        }
        else{
            loop.step();
        }
    }
}

void performIntergridTransfer(DA * oldDA, DA * newDA, DistTREE & oldDistTree, DistTREE & newDistTree, DistTREE & surrDistTree,double *& array, int ndof = 1){
  DA * surrDA = new DA(surrDistTree, MPI_COMM_WORLD, oldDA->getElementOrder());
  static std::vector<VECType> fineGhosted, surrGhosted;
  newDA->template createVector<VECType>(fineGhosted, false, true, ndof);
  surrDA->template createVector<VECType>(surrGhosted, false, true, ndof);
  std::fill(fineGhosted.begin(), fineGhosted.end(), 0);
  VECType *fineGhostedPtr = fineGhosted.data();
  VECType *surrGhostedPtr = surrGhosted.data();

  ot::distShiftNodes(*oldDA, array,
                     *surrDA, surrGhostedPtr + ndof * surrDA->getLocalNodeBegin(),
                     ndof);
  surrDA->template readFromGhostBegin<VECType>(surrGhostedPtr, ndof);
  surrDA->template readFromGhostEnd<VECType>(surrGhostedPtr, ndof);

  fem::MeshFreeInputContext<VECType, TREENODE>
    inctx{surrGhostedPtr,
          surrDA->getTNCoords(),
          (unsigned) surrDA->getTotalNodalSz(),
          &(*surrDistTree.getTreePartFiltered().cbegin()),
          surrDistTree.getTreePartFiltered().size(),
          *surrDA->getTreePartFront(),
          *surrDA->getTreePartBack()};

  fem::MeshFreeOutputContext<VECType, TREENODE>
    outctx{fineGhostedPtr,
           newDA->getTNCoords(),
           (unsigned) newDA->getTotalNodalSz(),
           &(*newDistTree.getTreePartFiltered().cbegin()),
           newDistTree.getTreePartFiltered().size(),
           *newDA->getTreePartFront(),
           *newDA->getTreePartBack()};

  // Hack for the current version
  std::vector<char> outDirty;
  newDA->template createVector<char>(outDirty, false, true, 1);
  newDA->setVectorByScalar(&(*outDirty.begin()), (std::array<char, 1>{0}).data(), false, true, 1);

  const RefElement *refel = newDA->getReferenceElement();
  fem::locIntergridTransfer(inctx, outctx, ndof, refel, &(*outDirty.begin()));

  newDA->template writeToGhostsBegin<VECType>(fineGhostedPtr, ndof, &(*outDirty.cbegin()));
  newDA->template writeToGhostsEnd<VECType>(fineGhostedPtr, ndof, false, &(*outDirty.cbegin()));

  delete[] array;

  newDA->createVector(array,false,false,ndof);
  newDA->template ghostedNodalToNodalVec<VECType>(fineGhostedPtr, array, true, ndof);

}

void checkIntergridTransfer(const double *array, ot::DA<DIM> * octDA,const ot::DistTree<unsigned int, DIM> &distTree, const unsigned int ndof){
    double *ghostedArray;
    octDA->nodalVecToGhostedNodal(array,ghostedArray, false,ndof);
    octDA->readFromGhostBegin(ghostedArray,ndof);
    octDA->readFromGhostEnd(ghostedArray,ndof);
    const size_t sz = octDA->getTotalNodalSz();
    auto partFront = octDA->getTreePartFront();
    auto partBack = octDA->getTreePartBack();
    const auto tnCoords = octDA->getTNCoords();
    const unsigned int nPe = octDA->getNumNodesPerElement();
    ot::MatvecBase<DIM, PetscScalar> treeloop(sz, ndof, octDA->getElementOrder(), tnCoords, ghostedArray, &(*distTree.getTreePartFiltered().cbegin()), distTree.getTreePartFiltered().size(), *partFront, *partBack);
    bool testPassed = true;
    while (!treeloop.isFinished())
    {
        if (treeloop.isPre() && treeloop.subtreeInfo().isLeaf())
        {
            const double * nodeCoordsFlat = treeloop.subtreeInfo().getNodeCoords();
            const PetscScalar * nodeValsFlat = treeloop.subtreeInfo().readNodeValsIn();
            for(int i = 0; i < nPe; i++){
                double correctValue = 0;
                for(int dim = 0; dim < DIM; dim++){
                    correctValue += nodeCoordsFlat[DIM*i+dim];
                }
                for(int dof = 0; dof < ndof; dof++) {
                    double interpolatedValue = nodeValsFlat[i * ndof + dof];
                    if (fabs(interpolatedValue - correctValue) > 1E-6) {
                        std::cout << "Value at (" << nodeCoordsFlat[DIM * i + 0] << " ," << nodeCoordsFlat[DIM * i + 1]
                                  << ") = " << interpolatedValue << "\n";
                        testPassed = false;
                    }
                }

            }
            treeloop.next();
        }
        else
            treeloop.step();
    }
    bool gtestPassed;
    MPI_Reduce(&testPassed,&gtestPassed,1,MPI_CXX_BOOL,MPI_LAND,0,MPI_COMM_WORLD);
    int rank; MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    if(!rank) {
      if (gtestPassed) {
        std::cout << GRN << "TEST passed" << NRM << "\n";
      } else {
        std::cout << RED << "TEST failed" << NRM << "\n";
      }
    }
    delete [] ghostedArray;
    MPI_Barrier(MPI_COMM_WORLD);
}
int main(int argc, char * argv[]){
  static const char *varname[]{"u"};

  using DENDRITE_UINT = unsigned  int;
    using TREENODE = ot::TreeNode<DENDRITE_UINT, DIM>;
    PetscInitialize(&argc, &argv, NULL, NULL);
    _InitializeHcurve(DIM);
    int eleOrder = 1;
    unsigned int ndof = 1;
    MPI_Comm comm = MPI_COMM_WORLD;
    ot::DistTree<unsigned int, DIM> oldDistTree;
    {
        std::vector<ot::TreeNode<unsigned int, DIM>> treePart;
        ot::createRegularOctree(treePart, 5, comm);
        oldDistTree = ot::DistTree<unsigned int, DIM>(treePart, comm);
    }
    ot::DA<DIM> *oldDA = new ot::DA<DIM>(oldDistTree, comm, eleOrder);
    /// Set Vector by a function
    double * coarseVec;
    oldDA->template createVector<VECType>(coarseVec,false,false,ndof);
    std::function<void(const double *, double *)> functionPointer = [&](const double *x, double *var) {
        double sum = 0.0;
        for (int d = 0; d < DIM; ++d)
            sum += x[d];
        var[0] = sum;
    };
    oldDA->setVectorByFunction(coarseVec,functionPointer,false,false,ndof);
    /// Refinement Flags
    std::vector<ot::OCT_FLAGS::Refine> octFlags(oldDistTree.getTreePartFiltered().size(),ot::OCT_FLAGS::Refine::OCT_NO_CHANGE);
    generateRefinementFlags(oldDA,octFlags,oldDistTree);
    ot::DistTree<unsigned int, DIM> newDistTree;
    ot::DistTree<unsigned int, DIM> surrDistTree;
    {
        std::vector<ot::TreeNode<DENDRITE_UINT, DIM>> newTree;
        std::vector<ot::TreeNode<DENDRITE_UINT, DIM>> surrTree;
        ot::SFC_Tree<DENDRITE_UINT , DIM>::distRemeshWholeDomain(oldDistTree.getTreePartFiltered(), octFlags, newTree, 0.3, comm);
        surrTree = ot::SFC_Tree<DENDRITE_UINT , DIM>::getSurrogateGrid(ot::RemeshPartition::SurrogateInByOut, oldDistTree.getTreePartFiltered(), newTree, comm);
        newDistTree = ot::DistTree<unsigned int, DIM>(newTree, comm);
        surrDistTree = ot::DistTree<unsigned int, DIM>(surrTree, comm,ot::DistTree<unsigned int,DIM>::NoCoalesce);
    }
    ot::DA<DIM> *newDA = new ot::DA<DIM>(newDistTree, comm, eleOrder);
    std::cout << "Number of elements in OldDA " << oldDA->getLocalElementSz() << "\n";
    std::cout << "Number of elements in NewDA " << newDA->getLocalElementSz() << "\n";
    performIntergridTransfer(oldDA,newDA,oldDistTree,newDistTree,surrDistTree,coarseVec,1);
    checkIntergridTransfer(coarseVec,newDA,newDistTree,1);
    std::swap(newDA,oldDA);
    delete newDA;
    PetscFinalize();
}
