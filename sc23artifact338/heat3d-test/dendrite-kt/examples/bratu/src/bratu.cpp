//
// Created by maksbh on 9/18/19.
//

#include <iostream>
#include <DendriteUtils.h>
#include <point.h>
#include <TalyEquation.h>
#include <BTNodeData.h>
#include <BTEquation.h>
#include <BTInputData.h>
#include <PETSc/Solver/LinearSolver.h>
#include <PETSc/PetscUtils.h>
#include <IO/VTU.h>
#include <PETSc/IO/petscVTU.h>
#include <Traversal/Analytic.h>

using namespace PETSc;
int main(int argc, char *argv[]) {

  dendrite_init(argc, argv);
  int rank = TALYFEMLIB::GetMPIRank();

  /// read parameters from config.txt
  BTInputData idata;
  /// use config.txt first
  std::ifstream configFile("config.txt");
  DENDRITE_UINT eleOrder = 0;
  DENDRITE_UINT level = 0;
  bool mfree = false;

  if (configFile.good()) {
    if (!idata.ReadFromFile()) {  /// read from file named "config.txt"
      throw std::runtime_error("[ERR] Error reading input data, check the config file!");
    }
    if (!idata.CheckInputData()) {
      throw std::runtime_error("[ERR] Problem with input data, check the config file!");
    }
    eleOrder = idata.basisFunction;
    level = idata.mesh_def.refine_lvl;
    mfree = idata.mfree;
  } else if (argc < 3) {
    TALYFEMLIB::PrintStatus("Usage: ", argv[0], " eleOrder level mfree");
    return -1;
  } else {
    eleOrder = std::atoi(argv[1]);
    level = std::atoi(argv[2]);
    mfree = std::atoi(argv[3]);
  }

  TALYFEMLIB::PrintInfo("Total number of processor = ", TALYFEMLIB::GetMPISize());
  TALYFEMLIB::PrintInfo("size of DendroInt ", sizeof(DendroIntL));
  TALYFEMLIB::PrintInfo("size of PetscInt ", sizeof(PetscInt));

  TALYFEMLIB::PrintStatus("eleOrder ", eleOrder);
  TALYFEMLIB::PrintStatus("Level ", level);
  TALYFEMLIB::PrintStatus("Mfree ", mfree);
  TALYFEMLIB::PrintStatus("DIM =  ", DIM);

  /// create DA
  std::vector<TREENODE> treePart;
  DA *octDA = createRegularDA(treePart, level, eleOrder);


  /// Problem size
  DENDRITE_UINT ndof = 1;
  DomainInfo physDomain;
  if (not(configFile.good())) {
    physDomain.min.fill(0.0);
    physDomain.max.fill(0.5);
  } else {
    physDomain.min = idata.mesh_def.min;
    physDomain.max = idata.mesh_def.max;
  }

  DomainExtents domain(physDomain);
  /// Equation and solver setup
  auto btEq = new TalyEquation<BTEquation, BTNodeData>(octDA,treePart, domain, ndof);
  NonlinearSolver *btSolver = setNonLinearSolver(btEq, octDA, ndof, mfree);
  if (configFile.good()) {
    /// apply solver parameter from config.txt
    idata.solverOptionsBT.apply_to_petsc_options("-bt_");
    {
      SNES m_snes = btSolver->snes();
      SNESSetOptionsPrefix(m_snes, "bt_");
      SNESSetFromOptions(m_snes);
    }
  }

  const auto analytic_sol = [&](TALYFEMLIB::ZEROPTV pos, int dof, double currentTime) {
    return -log(1 + cos(M_PI * (pos.x())));
  };

  /// For non linear solve : The vector you need and the corresponding dof
  btEq->setVectors({VecInfo(PLACEHOLDER_GUESS, ndof, 0)});

  /// BC setup
  btSolver->setDirichletBoundaryCondition([&](const TALYFEMLIB::ZEROPTV pos, unsigned int nodeID) -> Boundary {
    Boundary b;
    static constexpr double eps = 1e-14;
    if ((fabs(pos.x() - physDomain.min[0]) < eps) ||
        (fabs(pos.x() - physDomain.max[0]) < eps)) {

      b.addDirichlet(0, analytic_sol(pos, 0, 0));
    }
    return b;
  });
  const char *varname[]{"bt"};

  /// Solve
  btSolver->solve();

  /// Print Files
  petscVectopvtu(octDA, treePart,btSolver->getCurrentSolution(), "bratu", varname, domain, false, false, ndof);

  VecInfo v(btSolver->getCurrentSolution(), 1, 0);
  Analytic btAnalytic(octDA, treePart, v,  analytic_sol, domain, 0);
  btAnalytic.getL2error();
  /// Save vector file (only for regression test)
  if (idata.dump_vec) {
    petscDumpFilesforRegressionTest(octDA, btSolver->getCurrentSolution(), "solution_vec.vec");
  }
  dendrite_finalize(octDA);

}