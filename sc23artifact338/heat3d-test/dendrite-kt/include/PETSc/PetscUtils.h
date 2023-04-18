//
// Created by maksbh on 9/18/19.
//

#ifndef DENDRITEKT_PETSCUTILS_H
#define DENDRITEKT_PETSCUTILS_H

#include <PETSc/Solver/LinearSolver.h>
#include <PETSc/Solver/NonLinearSolver.h>
#include <PETSc/VecUtils.h>

namespace PETSc {
  template<typename Equation, typename NodeData>
  LinearSolver *setLinearSolver(TalyEquation<Equation, NodeData> *talyEq,
                                ot::DA<DIM> *da,
                                DENDRITE_UINT ndof = 1,
                                bool mfree = true,
                                const bool oneShotAssembly = false) {

    LinearSolver *solver = new LinearSolver(da, talyEq->getMat(), talyEq->getVec(), talyEq->getDomainMin(),
                                            talyEq->getDomainMax(), ndof, mfree, oneShotAssembly);
    if (da->isActive()) {
      solver->init();
    }
    return solver;

  }

  template<typename Equation, typename NodeData>
  NonlinearSolver *setNonLinearSolver(TalyEquation<Equation, NodeData> *talyEq,
                                      ot::DA<DIM> *da,
                                      DENDRITE_UINT ndof = 1,
                                      bool mfree = true,
                                      const bool oneShotAssembly = false) {

    NonlinearSolver *solver = new NonlinearSolver(da, talyEq->getMat(), talyEq->getVec(), talyEq->getDomainMin(),
                                                  talyEq->getDomainMax(), ndof, mfree, oneShotAssembly);
    if (da->isActive()) {
      solver->init();
    }
    return solver;

  }

  static void petscDumpFilesforRegressionTest(DA *octDA, const Vec &v, const std::string file_prefix = "vec") {
    PetscViewer viewer;
    int ierr = PetscViewerBinaryOpen(octDA->getCommActive(), file_prefix.c_str(), FILE_MODE_WRITE, &viewer);
    assert(ierr == 0);
    ierr = VecView(v, viewer);
    assert(ierr == 0);
    ierr = PetscViewerDestroy(&viewer);
    assert(ierr == 0);
  }


  static void recombineVec(DA *da, std::vector<VecInfo> &inVecs, VecInfo &outVec, bool isAllocated = false) {
    DENDRITE_UINT ndof = 0;
    for (auto &vec:inVecs) {
      ndof += vec.ndof;
      VecGetArrayRead(vec.v, &vec.val);
    }
    if(ndof != outVec.ndof){
      throw std::runtime_error("[Gather Error] : Dof did not match");
    }

    if (not(isAllocated)) {
      da->petscCreateVector(outVec.v, false, false, ndof);
    }
    PetscScalar *syncArray;
    VecGetArray(outVec.v, &syncArray);
    DENDRITE_UINT currentDof(0);
    const DENDRITE_UINT localNodalSz = da->getLocalNodalSz();
    for (DENDRITE_UINT node = 0; node < localNodalSz; node++) {
      currentDof = 0;
      for (const auto &vec:inVecs) {
        std::memcpy(&syncArray[node * ndof + currentDof], &vec.val[node * vec.ndof], sizeof(PetscScalar) * vec.ndof);
        currentDof += vec.ndof;
      }
    }

    for (auto &vec:inVecs) {
      VecRestoreArrayRead(vec.v, &vec.val);
    }
    VecRestoreArray(outVec.v, &syncArray);
  }

  static void recombineVec(DA *da, std::vector<VecUtils> &inVecs, VecUtils &outVec, bool isAllocated = false) {
    DENDRITE_UINT ndof = 0;
    for (auto &vec:inVecs) {
      ndof += vec.ndof;
      VecGetArrayRead(vec.petscVec, &vec.val);
    }

    if(ndof != outVec.ndof){
      throw std::runtime_error("[Gather Error] : Dof did not match");
    }

    if (not(isAllocated)) {
      da->petscCreateVector(outVec.petscVec, false, false, ndof);
    }
    PetscScalar *syncArray;
    VecGetArray(outVec.petscVec, &syncArray);
    DENDRITE_UINT currentDof(0);
    const DENDRITE_UINT localNodalSz = da->getLocalNodalSz();
    for (DENDRITE_UINT node = 0; node < localNodalSz; node++) {
      currentDof = 0;
      for (const auto &vec:inVecs) {
        std::memcpy(&syncArray[node * ndof + currentDof], &vec.val[node * vec.ndof], sizeof(PetscScalar) * vec.ndof);
        currentDof += vec.ndof;
      }
    }

    for (auto &vec:inVecs) {
      VecRestoreArrayRead(vec.petscVec, &vec.val);
    }
    VecRestoreArray(outVec.petscVec, &syncArray);
  }

  static void scatterVec(DA *octDA, std::vector<VecInfo> &outVec, const VecInfo &inVec, bool isAllocated = false) {
    DENDRITE_UINT ndof = 0;
    for (auto &vec:outVec) {
      ndof += vec.ndof;
    }
    if(ndof != inVec.ndof){
      throw std::runtime_error("[Scatter Error] Dof did not match");
    }
    if (not(isAllocated)) {
      for (auto &vec: outVec) {
        octDA->petscCreateVector(vec.v, false, false, vec.ndof);
      }
    }



    if (octDA->isActive()) {
      const PetscScalar *array;
      std::vector<PetscScalar *>arrayWrite(outVec.size());
      for(int i = 0; i < outVec.size(); i++){
        VecGetArray(outVec[i].v,&arrayWrite[i]);
      }
      VecGetArrayRead(inVec.v, &array);
      const int localNodalSz = octDA->getLocalNodalSz();
      const int totalDof = inVec.ndof;
      for (int i = 0; i < localNodalSz; i++) {
        int currentDofIndex = 0;
        for (int j = 0; j < outVec.size(); j++) {
          std::memcpy(&arrayWrite[j][i * outVec[j].ndof], &array[i * totalDof + currentDofIndex],
                      sizeof(PetscScalar) * outVec[j].ndof);
          currentDofIndex += outVec[j].ndof;
        }
      }
      for(int i = 0; i < outVec.size(); i++){
        VecRestoreArray(outVec[i].v,&arrayWrite[i]);
      }
      VecRestoreArrayRead(inVec.v, &array);
    }
  }

  static void scatterVec(DA *octDA, std::vector<VecUtils> &outVec, const VecUtils &inVec, bool isAllocated = false) {
    DENDRITE_UINT ndof = 0;
    for (auto &vec:outVec) {
      ndof += vec.ndof;
    }
    if(ndof != inVec.ndof){
      throw std::runtime_error("[Scatter Error] Dof did not match");
    }
    if (not(isAllocated)) {
      for (auto &vec: outVec) {
        octDA->petscCreateVector(vec.petscVec, false, false, vec.ndof);
      }
    }
    const PetscScalar *array;
    PetscScalar *arrayWrite;
    if (octDA->isActive()) {
      VecGetArrayRead(inVec.petscVec, &array);
      const int localNodalSz = octDA->getLocalNodalSz();
      const int totalDof = inVec.ndof;
      for (int i = 0; i < localNodalSz; i++) {
        int currentDofIndex = 0;
        for (int j = 0; j < outVec.size(); j++) {
          VecGetArray(outVec[j].petscVec, &arrayWrite);
          std::memcpy(&arrayWrite[i * outVec[j].ndof], &array[i * totalDof + currentDofIndex],
                      sizeof(PetscScalar) * outVec[j].ndof);
          VecRestoreArray(outVec[j].petscVec, &arrayWrite);
          currentDofIndex += outVec[j].ndof;
        }
      }
      VecRestoreArrayRead(inVec.petscVec, &array);
    }
  }

/**
 * @brief constructs to update equation and solver in case of AMR
 * @tparam Equation equationClass
 * @tparam NodeData NodeData
 * @tparam Args args
 * @param newDA newDA created after refinement
 * @param treePart treePartition
 * @param talyEquation talyEquation class. Old will be deleted and new will be created
 * @param domainExtents domainExtents.
 * @param subDomainBoundary required for surface Assembly. Null if not needed
 * @param ndof
 * @param solver solver class
 * @param args
 */
  template<typename Equation, typename NodeData, typename ... Args>
  static void updateEquationAndSolver(DA *newDA, const std::vector<TREENODE> &treePart,
                                      TalyEquation<Equation, NodeData> *&talyEquation,
                                      const DomainExtents &domainExtents,
                                      LinearSolver *&solver,
                                      const TimeInfo *timeInfo = nullptr,
                                      SubDomainBoundary *subDomainBoundary = nullptr,
                                      const int ndof = 1,
                                      Args ... args) {
    bool m_surfAssembly = talyEquation->performSurfAssembly();
    delete talyEquation;

    talyEquation = new TalyEquation<Equation, NodeData>(newDA, treePart, domainExtents, ndof, timeInfo, m_surfAssembly,
                                                        subDomainBoundary, args...);
    bool mfree = solver->get_mfree();
    bool oneShotAssembly = solver->oneShotAssembly();
    const auto bnd = solver->getBoundaryCondition();
    delete solver;
    solver = new LinearSolver(newDA, talyEquation->getMat(), talyEquation->getVec(), talyEquation->getDomainMin(),
                              talyEquation->getDomainMax(), ndof, mfree, oneShotAssembly);
    solver->setDirichletBoundaryCondition(bnd);
    if (newDA->isActive()) {
      solver->init();
    }
  }


  template<typename Equation, typename NodeData, typename ... Args>
  static void updateEquationAndSolver(DA *newDA, const std::vector<TREENODE> &treePart,
                                      TalyEquation<Equation, NodeData> *&talyEquation,
                                      const DomainExtents &domainExtents,
                                      NonlinearSolver *&solver,
                                      const TimeInfo *timeInfo = nullptr,
                                      SubDomainBoundary *subDomainBoundary = nullptr,
                                      const int ndof = 1,
                                      Args ... args) {
    bool m_surfAssembly = talyEquation->performSurfAssembly();
    delete talyEquation;
    talyEquation = new TalyEquation<Equation, NodeData>(newDA, treePart, domainExtents, ndof, timeInfo, m_surfAssembly,
                                                        subDomainBoundary, args...);
    bool mfree = solver->get_mfree();
    bool oneShotAssembly = solver->oneShotAssembly();
    const auto bnd = solver->getBoundaryCondition();
    delete solver;
    solver = new NonlinearSolver(newDA, talyEquation->getMat(), talyEquation->getVec(), talyEquation->getDomainMin(),
                                 talyEquation->getDomainMax(), ndof, mfree, oneShotAssembly);
    solver->setDirichletBoundaryCondition(bnd);
    if (newDA->isActive()) {
      solver->init();
    }
  }

  template<typename Equation, typename NodeData, typename ... Args>
  static void updateEquationAndSolver(TalyMesh<NodeData> *mesh, DA *newDA, const std::vector<TREENODE> &treePart,
                                      TalyEquation<Equation, NodeData> *&talyEquation,
                                      const DomainExtents &domainExtents,
                                      LinearSolver *&solver,
                                      const TimeInfo *timeInfo = nullptr,
                                      SubDomainBoundary *subDomainBoundary = nullptr,
                                      const int ndof = 1,
                                      Args ... args) {
    bool m_surfAssembly = talyEquation->performSurfAssembly();
    delete talyEquation;
    talyEquation = new TalyEquation<Equation, NodeData>(mesh, newDA, treePart, domainExtents, ndof, timeInfo,
                                                        m_surfAssembly, subDomainBoundary, args...);
    bool mfree = solver->get_mfree();
    bool oneShotAssembly = solver->oneShotAssembly();
    const auto bnd = solver->getBoundaryCondition();
    delete solver;
    solver = new LinearSolver(newDA, talyEquation->getMat(), talyEquation->getVec(), talyEquation->getDomainMin(),
                              talyEquation->getDomainMax(), ndof, mfree, oneShotAssembly);
    solver->setDirichletBoundaryCondition(bnd);
    if (newDA->isActive()) {
      solver->init();
    }
  }


  template<typename Equation, typename NodeData, typename ... Args>
  static void updateEquationAndSolver(TalyMesh<NodeData> *mesh, DA *newDA, const std::vector<TREENODE> &treePart,
                                      TalyEquation<Equation, NodeData> *&talyEquation,
                                      const DomainExtents &domainExtents,
                                      NonlinearSolver *&solver,
                                      const TimeInfo *timeInfo = nullptr,
                                      SubDomainBoundary *subDomainBoundary = nullptr,
                                      const int ndof = 1,
                                      Args ... args) {
    bool m_surfAssembly = talyEquation->performSurfAssembly();
    delete talyEquation;
    talyEquation = new TalyEquation<Equation, NodeData>(mesh, newDA, treePart, domainExtents, ndof, timeInfo,
                                                        m_surfAssembly, subDomainBoundary, args...);
    bool mfree = solver->get_mfree();
    bool oneShotAssembly = solver->oneShotAssembly();

    const auto bnd = solver->getBoundaryCondition();

    delete solver;
    solver = new NonlinearSolver(newDA, talyEquation->getMat(), talyEquation->getVec(), talyEquation->getDomainMin(),
                                 talyEquation->getDomainMax(), ndof, mfree, oneShotAssembly);
    solver->setDirichletBoundaryCondition(bnd);
    if (newDA->isActive()) {
      solver->init();
    }
  }
}

#endif //DENDRITEKT_PETSCUTILS_H
