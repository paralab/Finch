//
// Created by maksbh on 9/18/19.
//

#include <iostream>
#include <DendriteUtils.h>
#include <point.h>
#include <TalyEquation.h>
#include <SSHTNodeData.h>
#include <SSHTEquation.h>
#include <SSHTInputData.h>
#include <PETSc/Solver/LinearSolver.h>
#include <PETSc/PetscUtils.h>
#include <IO/VTU.h>
#include <PETSc/IO/petscVTU.h>
#include <Traversal/Analytic.h>

using namespace PETSc;
int main(int argc, char *argv[]) {
  dendrite_init(argc, argv);
  /// read parameters from config.txt
  SSHTInputData idata;

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

  /// Create DA
  std::vector<TREENODE> treePart;
  DA *octDA = createRegularDA(treePart, level, eleOrder);

  const char *varname[]{"T"};

  static const DENDRITE_UINT ndof = SSHTNodeData::valueno();

  /// Physical dimensions.
  DomainInfo physDomain;
  if (not(configFile.good())) {
    physDomain.min.fill(0.0);
    physDomain.max.fill(1.0);
  } else {
    physDomain.min = idata.mesh_def.min;
    physDomain.max = idata.mesh_def.max;
  }

  DomainExtents domainExtents(physDomain);

  /// Equation setup
  auto sshtEq = new TalyEquation<SSHTEquation, SSHTNodeData>(octDA, treePart, domainExtents);
  LinearSolver *sshtSolver = setLinearSolver(sshtEq, octDA, ndof, mfree);
  if (configFile.good()) {
    /// apply solver parameter from config.txt
    idata.solverOptionsSSHT.apply_to_petsc_options("-ssht_");
    {
      KSP m_ksp = sshtSolver->ksp();
      KSPSetOptionsPrefix(m_ksp, "ssht_");
      KSPSetFromOptions(m_ksp);
    }
  }
  /// Boundary condition
  sshtSolver->setDirichletBoundaryCondition([&](const TALYFEMLIB::ZEROPTV pos, unsigned int nodeID) -> Boundary {
    Boundary b;
    static constexpr double eps = 1e-14;
    double x = pos.x();
    double y = pos.y();

    bool on_wall = (fabs(x - physDomain.min[0]) < eps) ||
        (fabs(x - physDomain.max[0]) < eps) ||
        (fabs(y - physDomain.min[1]) < eps) ||
        (fabs(y - physDomain.max[1]) < eps);
#if (DIM >= 3)
    double z = pos.z();
    on_wall = on_wall || (fabs(z - physDomain.min[2]) < eps) || (fabs(z - physDomain.max[2]) < eps);
#endif
#if (DIM == 4)
    double t = pos.t();
    on_wall = on_wall ||  (fabs(t) < eps) || (fabs(t - physDomain.min[3]) < eps);
#endif
    if (on_wall) {
      b.addDirichlet(0, 0);
    }
    return b;
  });

  /// Solve
  sshtSolver->solve();



  /// Print files
  petscVectopvtu(octDA, treePart, sshtSolver->getCurrentSolution(), "ssht", varname, domainExtents, false, false, ndof);

  /// L2 error
  const auto analytic_sol = [](const TALYFEMLIB::ZEROPTV &pos, const DENDRITE_UINT dof, const DENDRITE_REAL time) {
#if (DIM == 2)
    return sin(M_PI * pos.x()) * sin(M_PI * pos.y());
#elif (DIM == 3)
    return sin(M_PI * pos.x()) * sin(M_PI * pos.y()) * sin(M_PI * pos.z());
#endif
  };
  VecInfo v(sshtSolver->getCurrentSolution(), 1, 0);
  Analytic sshtAnalytic(octDA, treePart, v, analytic_sol, physDomain);
  sshtAnalytic.getL2error();

  /// Save vector file (only for regression test)
  if (idata.dump_vec) {
    petscDumpFilesforRegressionTest(octDA, sshtSolver->getCurrentSolution(), "solution_vec.vec");
  }

  delete sshtEq;
  delete sshtSolver;
  dendrite_finalize(octDA);

}