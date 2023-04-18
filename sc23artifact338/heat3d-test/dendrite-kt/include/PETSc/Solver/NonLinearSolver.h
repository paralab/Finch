//
// Created by maksbh on 5/19/20.
//

#ifndef DENDRITEKT_NONLINEARSOLVER_H
#define DENDRITEKT_NONLINEARSOLVER_H
#include <sstream>

#include "Solver.h"
namespace PETSc {
class NonlinearSolver : public Solver {
  bool isAssembled = false;
 public:
  NonlinearSolver(DA* da,feMat<DIM> * mat, feVec<DIM> * vec, const Point<DIM> & domainMin,const Point<DIM> & domainMax, const DENDRITE_UINT dof, const bool mfree, const bool oneShotMatrixAssembly);

  ~NonlinearSolver() override;
  PetscErrorCode init() override;
  PetscErrorCode solve() override;
  virtual PetscErrorCode cleanup();

  PetscErrorCode jacobianMatMult(Vec in, Vec out) override;


  inline SNES & snes() {
    return m_snes;
  }
  inline bool get_mfree(){
    return m_matrixFree;
  }
  inline int getDof(){
    return m_uiDof;
  }

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
  static PetscErrorCode FormJacobian(SNES snes, Vec sol, Mat jac, Mat precond_matrix, void *ctx);
  static PetscErrorCode FormFunction(SNES snes, Vec in, Vec out, void *ctx);


  SNES m_snes = NULL;
  Vec m_guess = NULL;  // last guess passed into FormJacobian

};

}

#endif //DENDRITEKT_NONLINEARSOLVER_H
