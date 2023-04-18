//
// Created by maksbh on 11/8/20.
//
#include <iostream>
#include <distTree.h>
#include <oda.h>
#include <point.h>
#include <sfcTreeLoop_matvec_io.h>
#include <octUtils.h>

#include <intergridTransfer.h>

static constexpr unsigned int DIM = 2;
typedef ot::TreeNode<unsigned int, DIM> TREENODE;
typedef ot::DistTree<unsigned int, DIM> DistTREE;
constexpr unsigned int nchild = 1u << DIM;
/**
 * @brief Generate Flags that refine all except the boundary elements.
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
            if(not(boundaryOctant)){
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
    if(testPassed){
        std::cout << GRN << "TEST passed" << NRM << "\n";
    }
    else{
        std::cout << RED << "TEST failed" << NRM << "\n";
    }
}

ibm::Partition DomainDecider1(const double *physCoords, double physSize) {

    static bool XMortonOrder[nchild]{false, true, false, true};
    static bool YMortonOrder[nchild]{false, false, true, true};


    double coords[DIM];
    std::array<bool, nchild> isOutsideDomain; // isOutside = true means that we want to retain the points


    isOutsideDomain.fill(true);

    // Check node by node basis
    for (int n = 0; n < nchild; n++) {
        coords[0] = (physCoords[0] + XMortonOrder[n] * physSize) ;
        coords[1] = (physCoords[1] + YMortonOrder[n] * physSize) ;
        if((coords[0] >= 0.5) or (coords[0] == 0) or (coords[1] == 1.0) or (coords[1] == 0.0)){
            isOutsideDomain[n] = false;
        }

    }

    unsigned int numOutsidePoints = std::accumulate(isOutsideDomain.begin(), isOutsideDomain.end(), 0);

    if (numOutsidePoints == 0) { // No points is outside
        return ibm::IN;
    } else if (numOutsidePoints == nchild) { // All points are outside
        return ibm::OUT;
    }
    return ibm::INTERCEPTED;
}

int main(int argc, char *argv[]) {
    typedef unsigned int DENDRITE_UINT;
    PetscInitialize(&argc, &argv, NULL, NULL);
    _InitializeHcurve(DIM);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    m_uiMaxDepth = 25;
    const DENDRITE_UINT eleOrder = 1;
    if (argc < 2) {
        if (not(rank)) {
            std::cout << "Usage: level \n";
        }
        exit(EXIT_FAILURE);
    }

    const DENDRITE_UINT level = static_cast<DENDRITE_UINT>(std::atoi(argv[1]));

    MPI_Comm comm = MPI_COMM_WORLD;


    using DTree = ot::DistTree<unsigned int, DIM>;
    DTree oldDistTree = DTree::constructSubdomainDistTree(level, DomainDecider1,
                                                       comm);
    const std::vector<ot::TreeNode<unsigned int, DIM>> &treePart = oldDistTree.getTreePartFiltered();


    ot::DA<DIM> *oldDA = new ot::DA<DIM>(oldDistTree, comm, eleOrder);
    std::cout << oldDA->getLocalElementSz() << "\n";

    std::vector<VECType> coarseVec;
    const int ndof = 1;
    oldDA->template createVector<VECType>(coarseVec,false,false,ndof);
    std::function<void(const double *, double *)> functionPointer = [&](const double *x, double *var) {
        double sum = 0.0;
        for (int d = 0; d < DIM; ++d)
            sum += x[d];
        var[0] = sum;
    };
    oldDA->setVectorByFunction(coarseVec.data(),functionPointer,false,false,ndof);
    static const char * varname[]{"u"};


    /// Refinement Flags
    std::vector<ot::OCT_FLAGS::Refine> octFlags(oldDistTree.getTreePartFiltered().size(),ot::OCT_FLAGS::Refine::OCT_NO_CHANGE);
    generateRefinementFlags(oldDA,octFlags,oldDistTree);


    /// Intergrid transfer
    ot::DistTree<unsigned int, DIM> newDistTree;
    ot::DistTree<unsigned int, DIM> surrDistTree;
    {
        std::vector<ot::TreeNode<DENDRITE_UINT, DIM>> newTree;
        std::vector<ot::TreeNode<DENDRITE_UINT, DIM>> surrTree;
        DistTREE::distRemeshSubdomain(oldDistTree, octFlags, newDistTree, surrDistTree, 0.3);

    }
    ot::DA<DIM> *newDA = new ot::DA<DIM>(newDistTree, comm, eleOrder);
    ot::DA<DIM> *surrDA = new ot::DA<DIM>(surrDistTree, comm, eleOrder);
    std::cout << "Number of elements in OldDA " << oldDA->getLocalElementSz() << "\n";
    std::cout << "Number of elements in NewDA " << newDA->getLocalElementSz() << "\n";

    /// Intergrid Transfer
    static std::vector<VECType> fineGhosted, surrGhosted;
    newDA->template createVector<VECType>(fineGhosted, false, true, ndof);
    surrDA->template createVector<VECType>(surrGhosted,false, true, ndof);
    std::fill(fineGhosted.begin(), fineGhosted.end(), 0);
    VECType *fineGhostedPtr = fineGhosted.data();
    VECType *surrGhostedPtr = surrGhosted.data();

    // 1. Copy input data to ghosted buffer.
    ot::distShiftNodes(*oldDA,   coarseVec.data(),
                       *surrDA,     surrGhostedPtr + ndof * surrDA->getLocalNodeBegin(),
                       ndof);
    surrDA->template readFromGhostBegin<VECType>(surrGhostedPtr, ndof);
    surrDA->template readFromGhostEnd<VECType>(surrGhostedPtr, ndof);
    fem::MeshFreeInputContext<VECType, TREENODE>
        inctx{ surrGhostedPtr,
               surrDA->getTNCoords(),
               (unsigned) surrDA->getTotalNodalSz(),
               &(*surrDistTree.getTreePartFiltered().cbegin()),
               surrDistTree.getTreePartFiltered().size(),
               *surrDA->getTreePartFront(),
               *surrDA->getTreePartBack() };
    fem::MeshFreeOutputContext<VECType, TREENODE>
        outctx{fineGhostedPtr,
               newDA->getTNCoords(),
               (unsigned) newDA->getTotalNodalSz(),
               &(*newDistTree.getTreePartFiltered().cbegin()),
               newDistTree.getTreePartFiltered().size(),
               *newDA->getTreePartFront(),
               *newDA->getTreePartBack() };
    const RefElement * refel = newDA->getReferenceElement();
    fem::locIntergridTransfer(inctx, outctx, ndof, refel);
    newDA->template writeToGhostsBegin<VECType>(fineGhostedPtr, ndof);
    newDA->template writeToGhostsEnd<VECType>(fineGhostedPtr, ndof);
    double *newDAVec;
    newDA->createVector(newDAVec,false,false,ndof);
    newDA->template ghostedNodalToNodalVec<VECType>(fineGhostedPtr, newDAVec, true, ndof);
    checkIntergridTransfer(newDAVec,newDA,newDistTree,ndof);
    /// Bunch of stuff to delete
    PetscFinalize();

}