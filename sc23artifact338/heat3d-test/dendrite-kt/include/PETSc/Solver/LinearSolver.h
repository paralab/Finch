//
// Created by maksbh on 9/18/19.
//

#ifndef DENDRITEKT_LINEARSOLVER_H
#define DENDRITEKT_LINEARSOLVER_H

#include "Solver.h"

namespace PETSc {

class LinearSolver : public Solver {

  bool isMatAssembled = false;

 public:

  LinearSolver(DA* da,feMat<DIM> * mat, feVec<DIM> * vec, const Point<DIM> & domainMin,const Point<DIM> & domainMax, const DENDRITE_UINT dof, const bool mfree, const bool oneShotMatrixAssembly);
  /**
   * @brief initialize the solver object
   * @return error code.
   **/
  PetscErrorCode init() override;

  /**
   * @brief solve
   * @return error code
   */
  PetscErrorCode solve() override;

  /**
   * @brief performs the cleanup
   * @return error code
   */
  PetscErrorCode cleanup();

  /**
   * @brief matrix-free shell operations
   * @param [in] In constant Vec In
   * @param [out] Out The resultant vector out
   * @return error code.
   */

  PetscErrorCode jacobianMatMult(Vec In, Vec Out);


  /**
   *
   * @return the residual Vec
   */
  inline Vec VecGetResidual() {
    return m_vecRHS;
  }

  /**
   *
   * @return KSP cotext
   */
  inline KSP & ksp() {
    return m_ksp;
  }

  /**
   *
   * @return the dof for the solver.
   */
  inline int getDof() {
    return m_uiDof;
  }

  /**
   *
   * @return true if solver is matrix free, otherwise false
   */
  bool get_mfree();

  /**
   * @brief Destructor.
   */
  virtual ~LinearSolver();

 protected:
  using Solver::m_Mat;
  using Solver::m_Vec;
  using Solver::ShellMatMult;
  using Solver::m_octDA;
  using Solver::m_uiDof;
  using Solver::m_matrixFree;
  using Solver::m_matCompactRows;
#ifdef IBM
  using Solver::m_ibmDirichletNodes;
#endif
  KSP m_ksp = NULL;
};
}


#endif //DENDRITEKT_LINEARSOLVER_H
