//
// Created by maksbh on 9/18/19.
//

#ifndef DENDRITEKT_SOLVER_H
#define DENDRITEKT_SOLVER_H

#include <TalyEquation.h>
#include <TalyMat.h>
#include <TalyVec.h>

#include <vector>
#include <string>
#include <feMat.h>
#include <PETSc/BoundaryConditions.h>
#include <OctToPhysical.h>
#ifdef PROFILING
#include <Profiling/GlobalProfiler.h>
#endif

namespace PETSc{
class Solver {

 protected:
#ifdef IBM
  std::vector<PetscInt>  m_ibmDirichletNodes;
#endif

  ot::MatCompactRows m_matCompactRows;
  ot::DA<DIM> *m_octDA;

  feMat<DIM> * m_Mat = NULL;

  feVec<DIM> * m_Vec = NULL;


  Mat m_matJacobian;

  Vec m_vecRHS ;
  Vec m_vecSolution ;

  const bool m_matrixFree;


  const DENDRITE_UINT m_uiDof;
  BoundaryConditions boundary_conditions_;

  // for generating dirichlet bc
  std::function<Boundary(TALYFEMLIB::ZEROPTV, unsigned int)> m_boundaryCondition;



  OctToPhysical m_octToPhysical;

  const bool m_oneShotMatrixAssembly;

 public:

  Solver(DA* da,feMat<DIM> * mat, feVec<DIM> * vec,const Point<DIM> & domainMin,const Point<DIM> & domainMax, const DENDRITE_UINT dof,  const bool mfree, const bool oneShotMatrixAssembly);

  virtual ~Solver() {}

  /**
   * Initializes solver data structures (working vectors, SNES/KSP objects).
   * Should be called once after construction or a call to cleanup().
   * @return
   */
  virtual PetscErrorCode init() = 0;

  virtual PetscErrorCode solve() = 0;


  void setDirichletBoundaryCondition(const std::function<Boundary(const TALYFEMLIB::ZEROPTV, unsigned int)> &f);

  inline Vec  getCurrentSolution() const {
    if(m_octDA->isActive()) {
      assert(m_vecSolution != NULL);  // only valid after init()
    }
    return m_vecSolution;
  }

  inline BoundaryConditions& getBC() {
    return boundary_conditions_;
  }

  virtual PetscErrorCode jacobianMatMult(Vec _in, Vec _out) { throw std::runtime_error("Not implemented"); };



  /**
    *	@brief The Jacobian Matmult operation done a matrix free fashion
    *  @param _jac PETSC Matrix which is of shell type used in the time stepping
    *  @param _in  PETSC Vector which is the input vector
    *  @param _out PETSC Vector which is the output vector _jac*_in
    *  @return bool true if successful, false otherwise
    *
    *  See feMatrix.h for similar implementation
    **/

  static PetscErrorCode ShellMatMult(Mat M, Vec In, Vec Out) {
    Solver *contxt;
    MatShellGetContext(M, &contxt);
    return contxt->jacobianMatMult(In, Out);
  }

  virtual void updateBoundaries(ot::DA<DIM> *da);

  const std::function<Boundary(const TALYFEMLIB::ZEROPTV, unsigned int)> & getBoundaryCondition();

  inline bool oneShotAssembly()const{
    return m_oneShotMatrixAssembly;
  }
#ifdef IBM
  inline void setIBMDirichletNodes(std::vector<PetscInt> & dirichletNodes){
    m_ibmDirichletNodes.clear();
    m_ibmDirichletNodes.resize(dirichletNodes.size()*m_uiDof);
    for(int i = 0; i < dirichletNodes.size(); i++) {
      for (int dof = 0; dof < m_uiDof; dof++) {
        m_ibmDirichletNodes[i*m_uiDof+ dof] = dirichletNodes[i]*m_uiDof + dof;
      }
    }
    boundary_conditions_.assignIBMDirichletNodes(&m_ibmDirichletNodes);
  }
#endif
};
}
#endif //DENDRITEKT_SOLVER_H
