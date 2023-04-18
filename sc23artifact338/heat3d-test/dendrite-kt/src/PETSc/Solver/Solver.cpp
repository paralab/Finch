//
// Created by maksbh on 9/18/19.
//

#include <PETSc/Solver/Solver.h>
namespace PETSc {
Solver::Solver(DA* da,feMat<DIM> * mat, feVec<DIM> * vec, const Point<DIM> & domainMin,const Point<DIM> & domainMax,const DENDRITE_UINT dof,const bool mfree,const bool oneShotMatrixAssembly):
m_octDA(da),m_Mat(mat),m_Vec(vec),m_oneShotMatrixAssembly(oneShotMatrixAssembly),m_uiDof(dof),m_matrixFree(mfree),m_matCompactRows(da->getNumNodesPerElement(),dof,5000){


  // Matrix
  m_matJacobian = NULL;

  m_vecRHS = NULL;
  m_vecSolution = NULL;

  boundary_conditions_.setDof(dof);
  m_octToPhysical.init(domainMin,domainMax);
  boundary_conditions_.init(m_octToPhysical);

}

void Solver::setDirichletBoundaryCondition(const std::function<Boundary(const TALYFEMLIB::ZEROPTV, unsigned int)> &f) {
  m_boundaryCondition = f;
}

const std::function<Boundary(const TALYFEMLIB::ZEROPTV, unsigned int)> & Solver::getBoundaryCondition() {
  return m_boundaryCondition;
}


void Solver::updateBoundaries(ot::DA<DIM> *da) {
  if (m_boundaryCondition) {
    boundary_conditions_.clear();
    boundary_conditions_.addByNodalFunction(da, m_uiDof, m_boundaryCondition);
  }
}

}