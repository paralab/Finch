//
// Created by maksbh on 5/20/20.
//

#include <iostream>
#include <DendriteUtils.h>
#include <point.h>
#include <TalyEquation.h>
#include <HTEquation.h>
#include <HTNodeData.h>
#include <HTInputData.h>
#include <PETSc/Solver/LinearSolver.h>
#include <PETSc/PetscUtils.h>
#include <IO/VTU.h>
#include <PETSc/IO/petscVTU.h>
#include <Traversal/Analytic.h>
#include <unistd.h>
using namespace PETSc;
int main(int argc, char *argv[]) {

  dendrite_init(argc, argv);
  int rank = TALYFEMLIB::GetMPIRank();
  /// read parameters from config.txt
  HTInputData idata;

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

  std::vector<TREENODE> treePart;
  DA *octDA = createRegularDA(treePart, level, eleOrder);

  static const DENDRITE_UINT ndof = 1;
  const char *varname[]{"T"};
  /// Physical dimensions.
  DomainInfo physDomain;
  if (not(configFile.good())) {
    physDomain.min.fill(0.0);
    physDomain.max.fill(1.0);
  } else {
    physDomain.min = idata.mesh_def.min;
    physDomain.max = idata.mesh_def.max;
  }
  std::vector<DENDRITE_REAL> dt(1, 0.1);
  std::vector<DENDRITE_REAL> totalT(1, 1.0);
  TimeInfo ti(0.0, dt, totalT);
  DomainExtents domainExtents(physDomain);
  auto htEq = new TalyEquation<HTEquation, HTNodeData>(octDA, treePart,domainExtents, ndof,&ti, false,nullptr, &idata);

  LinearSolver *htSolver = setLinearSolver(htEq, octDA, ndof, mfree);
  if (configFile.good()) {
    /// apply solver parameter from config.txt
    idata.solverOptionsHT.apply_to_petsc_options("-ht_");
    {
      KSP m_ksp = htSolver->ksp();
      KSPSetOptionsPrefix(m_ksp, "ht_");
      KSPSetFromOptions(m_ksp);
    }
  }

  htSolver->setDirichletBoundaryCondition([&](const TALYFEMLIB::ZEROPTV pos, unsigned int nodeID) -> Boundary {
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
    if (on_wall) {
      b.addDirichlet(0, 0);
    }
    return b;
  });
  Vec prev_solution, prev_prev_solution;
  octDA->petscCreateVector(prev_solution, false, false, ndof);
  octDA->petscCreateVector(prev_prev_solution, false, false, ndof);
  std::function<void(const double *, double *)> initial_condition = [](const double *x, double *var) {
#if (DIM == 2)
    var[0] = sin(M_PI * x[0]) * sin(M_PI * x[1]);
#elif(DIM == 3)
    var[0] = sin(M_PI * x[0]) * sin(M_PI * x[1]) * sin(M_PI * x[2]);
#endif
  };
  octDA->petscSetVectorByFunction(prev_prev_solution, initial_condition, false, false, ndof);
  VecCopy(prev_prev_solution, prev_solution);
  htEq->setVectors({VecInfo(prev_solution, ndof, U_PRE),
                    VecInfo(prev_prev_solution, ndof, U_PRE_PRE)},
                   SYNC_TYPE::ALL);


  if (configFile.good()) {
    dt = idata.dt;
    totalT = idata.totalT;
  }
  TALYFEMLIB::PrintStatus("Pid = ", getpid());
  TALYFEMLIB::PrintStatus("Starting solve");

  while (ti.getCurrentTime() < ti.getEndTime() - 1e-15) {

    TALYFEMLIB::PrintStatus("Time = ", ti.getCurrentTime());
    htSolver->solve();
    VecCopy(prev_solution, prev_prev_solution);
    VecCopy(htSolver->getCurrentSolution(), prev_solution);
    ti.increment();
    /// Save solution when required
    if ((ti.getCurrentTime() >= (idata.OutputStartTime - 1e-16))
        && (ti.getTimeStepNumber() % idata.OutputInterval == 0)) {
      char fname[PATH_MAX];
      snprintf(fname, sizeof(fname), "%s_%03d", "ht", ti.getTimeStepNumber());
      petscVectopvtu(octDA, treePart,htSolver->getCurrentSolution(), fname, varname, domainExtents, false, false, ndof);
    }
  }

  const auto analytic_sol = [&](TALYFEMLIB::ZEROPTV pos, int dof, double currentTime) {
    if (configFile.good() && idata.forcing) {
#if(DIM == 2)
      return sin(M_PI * pos.x()) * sin(M_PI * pos.y()) * cos(M_PI * 2 * ti.getCurrentTime());
#elif (DIM == 3)
      return sin(M_PI * pos.x()) * sin(M_PI * pos.y()) * sin(M_PI * pos.z()) * cos(M_PI * 2 * ti.getCurrentTime());
#endif
    } else {
#if(DIM == 2)
      return sin(M_PI * pos.x()) * sin(M_PI * pos.y()) * exp(-currentTime);
#elif (DIM == 3)
      return sin(M_PI * pos.x()) * sin(M_PI * pos.y()) * sin(M_PI * pos.z()) * exp(-currentTime);
#endif
    }
  };
  /// Save vector file (only for regression test)
  if (idata.dump_vec) {
    petscDumpFilesforRegressionTest(octDA, htSolver->getCurrentSolution(), "solution_vec.vec");
  }

  VecInfo v(prev_solution, 1, 0);
  Analytic htAnalytic(octDA, treePart, v, analytic_sol, domainExtents, ti.getCurrentTime());
  htAnalytic.getL2error();

  /** Cleanup**/
  VecDestroy(&prev_solution);
  VecDestroy(&prev_prev_solution);
  delete htEq;
  delete htSolver;
  dendrite_finalize(octDA);

}