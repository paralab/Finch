//
// Created by maksbh on 11/28/18.
//

#include <PETSc/BoundaryConditions.h>
#include <talyfem/utils/utils.h>

namespace PETSc {

    void BoundaryConditions::clear() {
        m_boundaryRows.clear();
        m_boundaryValues.clear();
    }
    void BoundaryConditions::init(OctToPhysical &octToPhysical) {
        octToPhys_ = octToPhysical;
    }

    void BoundaryConditions::addByNodalFunction(ot::DA<DIM> *da, const  int ndof,
                       const std::function<Boundary(TALYFEMLIB::ZEROPTV, std::size_t)> &f) {

        if(!(da->isActive())){
            return;
        }

        std::vector<std::size_t > bdy_index;
        std::map<unsigned int, Boundary> bdy;



        da->getBoundaryNodeIndices(bdy_index);
        double coords[DIM];



        // Looping over the outer boundaries to check if the dirichlet BC needs to be specified
        const std::vector<DendroIntL> & localToGlobalMap = da->getNodeLocalToGlobalMap();
        for (int i = 0; i < bdy_index.size(); i++) {
          ot::treeNode2Physical(da->getTNCoords()[bdy_index[i] + da->getLocalNodeBegin()],da->getElementOrder(),coords);
          TALYFEMLIB::ZEROPTV zeroptvCoords;
            octToPhys_.convertCoordsToPhys(zeroptvCoords,coords);
            auto boundary = f(zeroptvCoords, bdy_index[i] + da->getLocalNodeBegin());
            bdy[bdy_index[i] + da->getLocalNodeBegin()] = boundary;
        }

        for (const auto &it : bdy) {
            for (const auto &direchlet : it.second.dirichlets) {
                PetscInt row = ndof * localToGlobalMap[it.first] + direchlet.first;

                m_boundaryRows.push_back(row);
                m_boundaryValues.push_back(direchlet.second);
            }

        }
    }



    PetscErrorCode BoundaryConditions::applyMatBC(ot::DA<DIM> *da, Mat mat) {

       /**
        * @author maksbh
        * Date July 29 2021:
        * Note for IBM: We assume that VecBC is always applied before the mat BC.
        * Boundary rows are updated in the applyResidualBC / apply Vec BC.
        * So, if we call boundaryRows.size() / boundaryRows.data(), this will contain IBM boundary values too.
        * No need to apply seperately.
        */

        int ierr;
        // TODO cache this if it works
        IS is;
        ierr = ISCreateGeneral(da->getCommActive(), m_boundaryRows.size(), m_boundaryRows.data(), PETSC_COPY_VALUES, &is);
        CHKERRQ(ierr);
        ierr = ISSortRemoveDups(is);CHKERRQ(ierr);
        ierr = MatZeroRowsIS(mat, is, 1.0, NULL, NULL); CHKERRQ(ierr);
        ISDestroy(&is);

        return 0;
    }

//    PetscErrorCode BoundaryConditions::applyMatIBMBC(ot::DA<DIM> *da, Mat mat, const std::vector<PetscInt> * dirichletNodeIDs ) {
//
//      int ierr;
//      // TODO cache this if it works
//      IS is;
//      ierr = ISCreateGeneral(da->getCommActive(), dirichletNodeIDs->size(), dirichletNodeIDs->data(), PETSC_COPY_VALUES, &is);
//      CHKERRQ(ierr);
//      ierr = ISSortRemoveDups(is);CHKERRQ(ierr);
//      ierr = MatZeroRowsIS(mat, is, 1.0, NULL, NULL); CHKERRQ(ierr);
//      ISDestroy(&is);
//
//      return 0;
//    }

    PetscErrorCode BoundaryConditions::applyVecBC(ot::DA<DIM>* da, Vec rhs) {

#ifdef IBM
        std::size_t oldBoundaryRowsSize = m_boundaryRows.size();
        m_boundaryRows.resize(oldBoundaryRowsSize + ibmDirichletNodes_->size());
        m_boundaryValues.resize(oldBoundaryRowsSize + ibmDirichletNodes_->size());
        std::memcpy(&m_boundaryRows[oldBoundaryRowsSize],ibmDirichletNodes_->data(), sizeof(PetscInt)*ibmDirichletNodes_->size());
        std::fill(m_boundaryValues.begin()+oldBoundaryRowsSize,m_boundaryValues.end(),0.0);
#endif
        int ierr = VecSetValues(rhs, m_boundaryRows.size(), m_boundaryRows.data(), m_boundaryValues.data(), INSERT_VALUES);
        VecAssemblyBegin(rhs);
        VecAssemblyEnd(rhs);
        return ierr;
    }

//    PetscErrorCode BoundaryConditions::applyVecIBMBC(ot::DA<DIM>* da, Vec rhs, const std::vector<PetscInt> & dirichletNodeIDs ) {
//
//      std::vector<PetscScalar> zeros(dirichletNodeIDs.size(),0.0);
//      int ierr = VecSetValues(rhs, dirichletNodeIDs.size(), dirichletNodeIDs.data(), zeros.data(), INSERT_VALUES);
//      VecAssemblyBegin(rhs);
//      VecAssemblyEnd(rhs);
//
//      return ierr;
//    }

    PetscErrorCode BoundaryConditions::applyResidualBC(ot::DA<DIM> *da, Vec residual) {

#ifdef IBM
      std::size_t oldBoundaryRowsSize = m_boundaryRows.size();
      m_boundaryRows.resize(oldBoundaryRowsSize + ibmDirichletNodes_->size());
      std::memcpy(&m_boundaryRows[oldBoundaryRowsSize],ibmDirichletNodes_->data(), sizeof(PetscInt)*ibmDirichletNodes_->size());
#endif
        std::vector<double> zeros(m_boundaryRows.size(), 0.0);
        int ierr = VecSetValues(residual, m_boundaryRows.size(), m_boundaryRows.data(), zeros.data(), INSERT_VALUES);
        VecAssemblyBegin(residual);
        VecAssemblyEnd(residual);

        return ierr;
    }

    PetscErrorCode BoundaryConditions::applyMatrixFreeBC(ot::DA<DIM> *da, Vec in, Vec out) {
        int ierr;
        IS is;

        ierr = ISCreateGeneral(da->getCommActive(), m_boundaryRows.size(), m_boundaryRows.data(), PETSC_USE_POINTER, &is);
        CHKERRQ(ierr);

        VecScatter scatter;
        ierr = VecScatterCreate(in, is, out, is, &scatter);
        CHKERRQ(ierr);

        ierr = VecScatterBegin(scatter, in, out, INSERT_VALUES, SCATTER_FORWARD);
        CHKERRQ(ierr);
        ierr = VecScatterEnd(scatter, in, out, INSERT_VALUES, SCATTER_FORWARD);
        CHKERRQ(ierr);

        ierr = VecScatterDestroy(&scatter);
        CHKERRQ(ierr);
        ierr = ISDestroy(&is);
        CHKERRQ(ierr);

        return 0;
    }

#ifdef IBM
//    PetscErrorCode BoundaryConditions::fillWithDirichletID_IBM(DA* da,Mat mat){
//      PetscInt start, end;
//      MatGetOwnershipRange(mat,&start,&end);
//      Vec v;
//      VecCreateMPI(da->getCommActive(), end-start, PETSC_DECIDE, &v);
//      MatGetDiagonal(mat,v);
//
//      const PetscScalar * array;
//      VecGetArrayRead(v,&array);
//      PetscInt lsz;
//      VecGetLocalSize(v,&lsz);
//      int counter = 0;
//      ibmDirichletNodes_.resize(lsz);
//      for(PetscInt i = 0; i < lsz; i++){
//        if(FEQUALS(array[i],0.0)){
//          ibmDirichletNodes_[counter++] = start + i;
//        }
//      }
//      ibmDirichletNodes_.resize(counter);
//      ibmDirichletNodes_.shrink_to_fit();
//      isAllocatedibmNodes = true;
//      VecRestoreArrayRead(v,&array);
//      VecDestroy(&v);
//      return 0;
//    }

//    PetscErrorCode BoundaryConditions::applyIBMBoundaryCondition(DA* da, Mat mat){
//      /// Do we need to add the vec also???
//      if(not isAllocatedibmNodes){
//        fillWithDirichletID_IBM(da,mat);
//      }
//      for (const auto &row : ibmDirichletNodes_) {
//        for(int dof = 0; dof < ndof_; dof++) {
//          MatSetValue(mat, row * ndof_ + dof, row * ndof_ + dof, 1.0, INSERT_VALUES);
//        }
//      }
//      MatAssemblyBegin(mat,MAT_FINAL_ASSEMBLY);
//      MatAssemblyEnd(mat,MAT_FINAL_ASSEMBLY);
//
//      return 0;
//    }
//
//  PetscErrorCode BoundaryConditions::applyIBMBoundaryCondition(DA* da, Mat mat, const std::vector<PetscInt> * dirichletNodeIDs){
//    for (const auto &row : *dirichletNodeIDs) {
//      for(int dof = 0; dof < ndof_; dof++) {
//        MatSetValue(mat, row * ndof_ + dof, row * ndof_ + dof, 1.0, INSERT_VALUES);
//      }
//    }
//    MatAssemblyBegin(mat,MAT_FINAL_ASSEMBLY);
//    MatAssemblyEnd(mat,MAT_FINAL_ASSEMBLY);
//    return 0;
//  }
#endif


}