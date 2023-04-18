//
// Created by maksbh on 9/18/19.
//

#include <iostream>
#include <DendriteUtils.h>
#include <point.h>
#include <TalyEquation.h>
#include <PPNodeData.h>
#include <PPEquation.h>
#include <PETSc/Solver/LinearSolver.h>
#include <PETSc/PetscUtils.h>
#include <IO/VTU.h>
#include <PETSc/IO/petscVTU.h>
//#include <PostProcessing/postProcessing.h>

using namespace PETSc;
int main(int argc, char *argv[]) {

  dendrite_init(argc, argv);
  int rank = TALYFEMLIB::GetMPIRank();
  if (argc < 3) {
    if (not(rank))
      std::cout << "Usage: " << argv[0]
                << " eleOrder level mfree"
                << std::endl;
    return 0;
  }

  TALYFEMLIB::PrintInfo("Total number of processor = ", TALYFEMLIB::GetMPISize());
  TALYFEMLIB::PrintInfo("size of DendroInt ", sizeof(DendroIntL));
  TALYFEMLIB::PrintInfo("size of PetscInt ", sizeof(PetscInt));

  const DENDRITE_UINT eleOrder = static_cast<DENDRITE_UINT>(std::atoi(argv[1]));
  const DENDRITE_UINT level = static_cast<DENDRITE_UINT>(std::atoi(argv[2]));
  const bool mfree = static_cast<bool>(std::atoi(argv[3]));
  TALYFEMLIB::PrintStatus("eleOrder ", eleOrder);
  TALYFEMLIB::PrintStatus("Level ", level);
  TALYFEMLIB::PrintStatus("Mfree ", mfree);
  TALYFEMLIB::PrintStatus("DIM =  ", DIM);

  std::vector<TREENODE> treePart;
  DA *octDA = createRegularDA(treePart, level, eleOrder);

  DENDRITE_UINT ndof = 2;

  /// Physical dimensions.
  DomainInfo physDomain;
  physDomain.min.fill(0.0);
  physDomain.max.fill(1.0);

  DomainExtents domainExtents(physDomain);

  auto ppEq = new TalyEquation<PPEquation, PPNodeData>(octDA, treePart, domainExtents, ndof);
  LinearSolver *ppSolver = setLinearSolver(ppEq, octDA, ndof, mfree);

  ppSolver->setDirichletBoundaryCondition([&](const TALYFEMLIB::ZEROPTV pos, unsigned int nodeID) -> Boundary {
    Boundary b;
    return b;
  });

  std::function<void(const double *, double *)> f_initial = [](const double *x, double *var) {
    var[0] = 2.0;
    var[1] = -1.0;
  };
  Vec prev_solution;
  octDA->petscCreateVector(prev_solution, false, false, ndof);
  octDA->petscSetVectorByFunction(prev_solution, f_initial, false, false, ndof);
  const char *varName[]{"u", "v"};
  petscVectopvtu(octDA, treePart, prev_solution, "pp_init", varName, domainExtents, false, false, ndof);
  ppEq->setVectors({VecInfo(prev_solution, ndof, U_PRE)});
  std::vector<DENDRITE_REAL> dt(1, 0.1);
  std::vector<DENDRITE_REAL> totalT(1, 0.2);
  TimeInfo ti(0.0, dt, totalT);
  ppEq->setTime(&ti);
  while (ti.getCurrentTime() < ti.getEndTime()) {
    ti.increment();
    ti.print();
    ppSolver->solve();
    VecCopy(ppSolver->getCurrentSolution(), prev_solution);
  }

  VecDestroy(&prev_solution);

  dendrite_finalize(octDA);

}