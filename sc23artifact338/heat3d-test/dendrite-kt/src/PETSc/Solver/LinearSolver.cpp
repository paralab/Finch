//
// Created by maksbh on 9/18/19.
//


#include <PETSc/Solver/LinearSolver.h>
#include <petsc.h>
#include <DendriteUtils.h>

namespace PETSc {
  LinearSolver::LinearSolver(DA *da, feMat<DIM> *mat, feVec<DIM> *vec, const Point<DIM> &domainMin,
                             const Point<DIM> &domainMax, const DENDRITE_UINT dof, const bool mfree,
                             const bool oneShotMatrixAssembly)
    : Solver(da, mat, vec, domainMin, domainMax, dof, mfree, oneShotMatrixAssembly) {

  }

  PetscErrorCode LinearSolver::init() {
    int ierr;  // for petsc error codes
    ierr = m_octDA->petscCreateVector(m_vecSolution, false, false, m_uiDof);
    CHKERRQ(ierr);
    ierr = m_octDA->petscCreateVector(m_vecRHS, false, false, m_uiDof);
    CHKERRQ(ierr);

    if (m_matrixFree) {
      int matsize;
      ierr = VecGetLocalSize(m_vecSolution, &matsize);
      CHKERRQ(ierr);
      ierr = MatCreateShell(m_octDA->getCommActive(), matsize, matsize, PETSC_DETERMINE, PETSC_DETERMINE, this,
                            &m_matJacobian);
      CHKERRQ(ierr);
      ierr = MatShellSetOperation(m_matJacobian, MATOP_MULT, (void (*)()) (ShellMatMult));
      CHKERRQ(ierr);
    } else {
      ierr = m_octDA->createMatrix(m_matJacobian, MATAIJ, m_uiDof);
      CHKERRQ(ierr);
      // avoids new nonzero errors when boundary conditions are used
      ierr = MatSetOption(m_matJacobian, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
      CHKERRQ(ierr);
      ierr = MatSetOption(m_matJacobian, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
      CHKERRQ(ierr);
    }


    // create a KSP context
    ierr = KSPCreate(m_octDA->getCommActive(), &m_ksp);
    CHKERRQ(ierr);
    ierr = KSPSetOperators(m_ksp, m_matJacobian, m_matJacobian);
    CHKERRQ(ierr);
    ierr = KSPSetFromOptions(m_ksp);
    CHKERRQ(ierr);
    return 0;
  }

  PetscErrorCode LinearSolver::solve() {

    if (m_octDA->isActive()) {
      int ierr;
#ifdef PROFILING
      PetscLogEventBegin(Profiling::vecAssembly,0,0,0,0);
#endif
      // update Boundaries
      updateBoundaries(m_octDA);
      ierr = VecZeroEntries(m_vecRHS);
      CHKERRQ(ierr);
#ifdef PROFILING
      PetscLogEventBegin(Profiling::vecElementalAssembly,0,0,0,0);
#endif
      m_Vec->computeVec(m_vecSolution, m_vecRHS);
#ifdef PROFILING
      PetscLogEventEnd(Profiling::vecElementalAssembly,0,0,0,0);
#endif
#ifdef PROFILING
      PetscLogEventBegin(Profiling::vecBC,0,0,0,0);
#endif

      // apply vec BC
      ierr = getBC().applyVecBC(m_octDA, m_vecRHS);
#ifdef PROFILING
      PetscLogEventEnd(Profiling::vecBC,0,0,0,0);
#endif
      CHKERRQ(ierr);
      ierr = VecAssemblyBegin(m_vecRHS);
      CHKERRQ(ierr);
      ierr = VecAssemblyEnd(m_vecRHS);
      CHKERRQ(ierr);
#ifdef PROFILING
      PetscLogEventEnd(Profiling::vecAssembly,0,0,0,0);
#endif

#ifdef PROFILING
      PetscLogEventBegin(Profiling::matAssembly,0,0,0,0);
#endif
      if (!m_matrixFree and not(isMatAssembled)) {


        ierr = MatZeroEntries(m_matJacobian);
#ifdef PROFILING
        PetscLogEventBegin(Profiling::matElementalAssembly,0,0,0,0);
#endif
        m_Mat->getAssembledMatrix(&m_matJacobian, MATAIJ,m_matCompactRows);
#ifdef PROFILING
        PetscLogEventEnd(Profiling::matElementalAssembly,0,0,0,0);
#endif
        CHKERRQ(ierr);

#ifdef IBM
        ierr = MatAssemblyBegin(m_matJacobian, MAT_FINAL_ASSEMBLY);
        CHKERRQ(ierr);
        ierr = MatAssemblyEnd(m_matJacobian, MAT_FINAL_ASSEMBLY);
        CHKERRQ(ierr);
#else

        ierr = MatAssemblyBegin(m_matJacobian, MAT_FINAL_ASSEMBLY);
        CHKERRQ(ierr);
        ierr = MatAssemblyEnd(m_matJacobian, MAT_FINAL_ASSEMBLY);
        CHKERRQ(ierr);
#endif
#ifdef PROFILING
        PetscLogEventBegin(Profiling::matBC,0,0,0,0);
#endif
        ierr = getBC().applyMatBC(m_octDA, m_matJacobian);
        CHKERRQ(ierr);
#ifdef PROFILING
        PetscLogEventEnd(Profiling::matBC,0,0,0,0);
#endif
        if (m_oneShotMatrixAssembly) {
          // For one shot assembly, once the matrix is assembled, we can delete the
          // memory for the matrix allocations as this will never be used.
          if(m_matCompactRows.getRowIdxs().capacity() != 0) {
            ot::MatCompactRows _rows(m_octDA->getNumNodesPerElement(), m_uiDof,0);
            std::swap(_rows, m_matCompactRows);
          }
          isMatAssembled = true;
        }
      }

#ifdef PROFILING
      PetscLogEventEnd(Profiling::matAssembly,0,0,0,0);
#endif
      KSPSetInitialGuessNonzero(m_ksp, PETSC_TRUE);
      ierr = KSPSolve(m_ksp, m_vecRHS, m_vecSolution);


      CHKERRQ(ierr);
    }
    return 0;

  }

  PetscErrorCode LinearSolver::jacobianMatMult(Vec In, Vec Out) {

    assert(m_matrixFree);
    int ierr;

    ierr = VecZeroEntries(Out);
    CHKERRQ(ierr);

    m_Mat->matVec(In, Out);
    ierr = getBC().applyMatrixFreeBC(m_octDA, In, Out);
    CHKERRQ(ierr);

    return 0;
  }

  LinearSolver::~LinearSolver() {
    cleanup();
  }

  PetscErrorCode LinearSolver::cleanup() {
    int ierr;
    ierr = KSPDestroy(&m_ksp);
    CHKERRQ(ierr);
    ierr = MatDestroy(&m_matJacobian);
    CHKERRQ(ierr);
    ierr = VecDestroy(&m_vecSolution);
    CHKERRQ(ierr);
    ierr = VecDestroy(&m_vecRHS);
    CHKERRQ(ierr);
    return 0;
  }

  bool LinearSolver::get_mfree() {
    return (m_matrixFree);
  }

}