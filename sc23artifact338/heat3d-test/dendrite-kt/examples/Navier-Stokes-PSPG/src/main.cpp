//
// Created by maksbh on 5/21/20.
//


#include <iostream>
#include <DendriteUtils.h>
#include <point.h>
#include <TalyEquation.h>
#include <NSNodeData.h>
#include <NSEquation.h>
#include <PETSc/Solver/LinearSolver.h>
#include <PETSc/PetscUtils.h>
#include <IO/VTU.h>
#include <Traversal/Analytic.h>
#include <PETSc/IO/petscVTU.h>
#include <NSUtils.h>
#include <NSBoundaryConditions.h>

using namespace PETSc;

int main(int argc, char *argv[]) {

  dendrite_init(argc, argv);
  int rank = TALYFEMLIB::GetMPIRank();

  NSInputData inputData;
  if (!(inputData.ReadFromFile())) {
    if (!rank) {
      throw TALYFEMLIB::TALYException() << "Can't read the config file \n";
    }
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  if (!inputData.CheckInputData()) {
    throw std::runtime_error("[ERR] Problem with input data, check the config file!");
  }
  TALYFEMLIB::PrintInfo("Total number of processor = ", TALYFEMLIB::GetMPISize());
  TALYFEMLIB::PrintInfo("size of DendroInt ", sizeof(DendroIntL));
  TALYFEMLIB::PrintInfo("size of PetscInt ", sizeof(PetscInt));

  if (argc == 4) {
    inputData.basisFunction = static_cast<TALYFEMLIB::kBasisFunction>(std::atoi(argv[1]));
    inputData.meshDef.baseLevel = static_cast<DENDRITE_UINT>(std::atoi(argv[2]));
    inputData.mfree = static_cast<bool>(std::atoi(argv[3]));
  }

  TALYFEMLIB::PrintStatus("eleOrder ", inputData.basisFunction);
  TALYFEMLIB::PrintStatus("Level ", inputData.meshDef.baseLevel);
  TALYFEMLIB::PrintStatus("Mfree ", inputData.mfree);
  TALYFEMLIB::PrintStatus("DIM =  ", DIM);

  std::vector<TREENODE> treePart;
  DA *octDA = createRegularDA(treePart, inputData.meshDef.baseLevel, inputData.basisFunction);

  DomainInfo physDomain;
  physDomain.min = inputData.meshDef.min;
  physDomain.max = inputData.meshDef.max;
  DomainExtents domain(physDomain);
  performRefinement(octDA, domain, treePart, inputData);

  std::vector<DENDRITE_REAL> dt = inputData.dt;
  std::vector<DENDRITE_REAL> totalT = inputData.totalT;
  if (argc == 4) {
    dt = std::vector<DENDRITE_REAL>(1, 0.1);
    totalT = std::vector<DENDRITE_REAL>(1, 0.2);
  }
  TimeInfo ti(0.0, dt, totalT);

  NSBoundaryConditions nsBC(&inputData, &ti);

  auto nsEq = new TalyEquation<NSEquation, NSNodeData>(octDA,treePart,domain, NSNodeData::NS_DOF,&ti, false, nullptr, &inputData);
  NonlinearSolver *nsSolver = setNonLinearSolver(nsEq, octDA, NSNodeData::NS_DOF, inputData.mfree);
  inputData.solverOptionsNS.apply_to_petsc_options("-");
  SNES m_snes = nsSolver->snes();
  SNESSetFromOptions(m_snes);
  if (inputData.ifMMS) {
    nsBC.setAnalyticalFunction(AnalyticalSolution);
  }

  nsSolver->setDirichletBoundaryCondition([&](const TALYFEMLIB::ZEROPTV &pos, int nodeID) -> Boundary {
    Boundary b;
    nsBC.getMomentumBoundaryCondition(b, pos);
    return b;
  });

  Vec prev1Solution, prev2Solution;
  octDA->petscCreateVector(prev1Solution, false, false, NSNodeData::NS_DOF);
  octDA->petscCreateVector(prev2Solution, false, false, NSNodeData::NS_DOF);
  setInitialCondition(octDA, inputData, prev1Solution);

  nsEq->setVectors({VecInfo(PLACEHOLDER_GUESS, NSNodeData::NS_DOF, NSNodeData::VEL_X),
                         VecInfo(prev1Solution,NSNodeData::NS_DOF,NSNodeData::VEL_X_PRE1),
                         VecInfo(prev2Solution,NSNodeData::NS_DOF,NSNodeData::VEL_X_PRE2)
  });

  VecCopy(prev1Solution, nsSolver->getCurrentSolution());
#if(DIM == 3)
  static const char *varname[]{"u", "v", "w", "p"};
#endif
#if(DIM == 2)
  static const char *varname[]{"u", "v", "p"};
#endif
  petscVectopvtu(octDA,treePart, nsSolver->getCurrentSolution(), "ns_init", varname, domain, false, false, NSNodeData::NS_DOF);
  double error[NSNodeData::NUM_VARS];
  while (ti.getCurrentTime() < ti.getEndTime() - 1e-15) {
    ti.increment();
    ti.print();
    nsSolver->solve();
    VecCopy(prev1Solution, prev2Solution);
    VecCopy(nsSolver->getCurrentSolution(), prev1Solution);
    if (inputData.ifMMS) {
      VecInfo v(prev1Solution, NSNodeData::NS_DOF, 0);
      Analytic NSAnalytic(octDA,treePart, v, AnalyticalSolution, physDomain, ti.getCurrentTime());
      NSAnalytic.getL2error(error);
      NSAnalytic.getL2error();
      if(not(rank)) {

        for (int i = 0; i < NSNodeData::NUM_VARS; i++) {
          std::cout << error[i] << " ";
        }
        std::cout << "\n";
      }
    }
  }

  petscVectopvtu(octDA,treePart, nsSolver->getCurrentSolution(), "ns", varname, domain, false, false, NSNodeData::NS_DOF);
  /// Save vector file (only for regression test)
  if (inputData.dump_vec) {
    petscDumpFilesforRegressionTest(octDA, nsSolver->getCurrentSolution(), "solution_vec.vec");

  }

  /** Clean up memory **/
  VecDestroy(&prev1Solution);
  VecDestroy(&prev2Solution);
  delete nsEq;
  delete nsSolver;

  dendrite_finalize(octDA);

}