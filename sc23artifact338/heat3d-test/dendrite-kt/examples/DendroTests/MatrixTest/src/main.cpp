//
// Created by maksbh on 10/29/20.
//

#include <PETSc/IO/petscVTU.h>
#include "matrix.h"

/**
 * This function only refines the elements marked as boundary octants
 * @param octDA
 * @param refineFlags
 */
void generateRefinementFlags(const ot::DA<DIM> *octDA, const std::vector<TREENODE> &treePart,
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

    ot::MatvecBaseCoords<DIM> loop(sz, eleOrder, false, 0, tnCoords, &(*treePart.cbegin()), treePart.size(), *partFront,
                                   *partBack);
    while (!loop.isFinished()) {
        if (loop.isPre() && loop.subtreeInfo().isLeaf()) {
          int numInside = 0;
          const auto & coords = loop.subtreeInfo().getNodeCoords();
          for(int i = 0; i < 8; i++){
            double dist = (coords[3*i + 0] - 0.5)*(coords[3*i + 0] - 0.5)+
                          (coords[3*i + 1] - 0.5)*(coords[3*i + 1] - 0.5)+(coords[3*i + 2] - 0.5)*(coords[3*i + 2] - 0.5);
            if(dist < (0.25*0.25)){
              numInside++;
            }
          }
          if((numInside != 0) and (numInside != 8)){
            refineFlags[counter] =  ot::OCT_FLAGS::Refine::OCT_REFINE;
          }
            counter++;
            loop.next();

        } else {
            loop.step();

        }
    }

}

int main(int argc, char *argv[]) {

    PetscInitialize(&argc, &argv, NULL, NULL);
    _InitializeHcurve(DIM);
    int eleOrder = 1;

    std::vector<ot::TreeNode<DENDRITE_UINT, DIM>> newTree;
    std::vector<ot::TreeNode<DENDRITE_UINT, DIM>> surrTree;
    std::vector<ot::TreeNode<unsigned int, DIM>> treePart;
    ot::createRegularOctree(treePart, 3, MPI_COMM_WORLD);
    DomainInfo domainInfo;
    domainInfo.min.fill(0.0);
    domainInfo.max.fill(1.0);
    DomainExtents domain(domainInfo);
    ot::DA<DIM> *newDA = new ot::DA<DIM>(treePart, MPI_COMM_WORLD, eleOrder);
    {
        std::vector<ot::OCT_FLAGS::Refine> octFlags;
        generateRefinementFlags(newDA, treePart, octFlags);
        ot::SFC_Tree<DENDRITE_UINT, DIM>::distRemeshWholeDomain(treePart, octFlags, newTree, surrTree, 0.3,
                                                                MPI_COMM_WORLD);
        DA *octDA = new DA(newTree, MPI_COMM_WORLD, eleOrder, 100, 0.3);
        std::swap(octDA, newDA);
        std::swap(newTree,treePart);
        delete octDA;
    }


    Vec funcVec;
    newDA->petscCreateVector(funcVec, false, false, 1);
    std::function<void(const double *, double *)> functionPointer = [&](const double *x, double *var) {
        var[0] = sin(M_PI * x[0]) * sin(M_PI * x[1]) * sin(M_PI * x[2]);
    };
    newDA->petscSetVectorByFunction(funcVec, functionPointer, false, false, 1);
    static const char *varname[]{"u"};
    PETSc::petscVectopvtu(newDA, treePart,funcVec, "kTVec", varname, domain, false, false, 1);
//

    matrix<DIM> mat(newDA, treePart, 0);
    Mat J;
    newDA->createMatrix(J, MATAIJ);
    MatSetOption(J, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    mat.getAssembledMatrix(&J, MATAIJ);
    MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY);
    Vec result, matVecResult;
    newDA->petscCreateVector(result, false, false, 1);
    newDA->petscCreateVector(matVecResult, false, false, 1);
    MatMult(J, funcVec, result);
    const double *value, *value1;
    VecGetArrayRead(result, &value);

    matrix<DIM> matVec(newDA, treePart, 0);
    VecSet(matVecResult, 0.0);
    matVec.matVec(funcVec, matVecResult);
    VecGetArrayRead(matVecResult, &value1);
    std::vector<double> maxDiff(newDA->getLocalNodalSz(),0.0);
    for (int i = 0; i < newDA->getLocalNodalSz(); i++) {
        maxDiff[i] = fabs(value[i] - value1[i]);
    }

    std::cout << *std::max_element(maxDiff.begin(),maxDiff.end()) << "\n";
    VecRestoreArrayRead(matVecResult, &value1);
    VecRestoreArrayRead(result, &value);
//
    PetscFinalize();
}