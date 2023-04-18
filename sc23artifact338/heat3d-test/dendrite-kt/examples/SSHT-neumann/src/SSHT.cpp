//
// Created by maksbh on 7/28/20.
//

//
// Created by maksbh on 9/18/19.
//

#include <iostream>
#include <DendriteUtils.h>
#include <point.h>
#include <TalyEquation.h>
#include <SSHTNodeData.h>
#include <SSHTEquation.h>
#include <PETSc/Solver/LinearSolver.h>
#include <PETSc/PetscUtils.h>
#include <IO/VTU.h>
#include <PETSc/IO/petscVTU.h>
#include <Traversal/Analytic.h>

using namespace PETSc;
int main(int argc, char *argv[]) {
  dendrite_init(argc, argv);
  if (argc < 3) {
    TALYFEMLIB::PrintStatus("Usage: ", argv[0], " eleOrder level mfree");
    return EXIT_FAILURE;
  }
  DENDRITE_UINT eleOrder = static_cast<DENDRITE_UINT>(std::atoi(argv[1]));
  DENDRITE_UINT  level   = static_cast<DENDRITE_UINT>(std::atoi(argv[2]));
  bool mfree = static_cast<bool>(std::atoi(argv[3]));


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

  static const DENDRITE_UINT ndof = 1;

  /// Physical dimensions.
  DomainInfo physDomain;
  physDomain.min.fill(0.0);
  physDomain.max.fill(1.0);

  DomainExtents domainExtents(physDomain);

  DomainBoundary domainBoundary(octDA,domainExtents);

  /// Equation setup : Note do not use dynamic_cast
  auto sshtEq = new TalyEquation<SSHTEquation, SSHTNodeData>(octDA,treePart, domainExtents,ndof,static_cast<SubDomainBoundary*>(&domainBoundary));
  LinearSolver *sshtSolver = setLinearSolver(sshtEq, octDA, ndof, mfree);
  /// Boundary condition
  sshtSolver->setDirichletBoundaryCondition([&](const TALYFEMLIB::ZEROPTV pos, unsigned int nodeID) -> Boundary {
    Boundary b;
    static constexpr double eps = 1e-14;
    double x = pos.x();
    double y = pos.y();

    bool on_wall = (fabs(x - physDomain.min[0]) < eps) or
                   (fabs(x - physDomain.max[0]) < eps);
    if (on_wall) {
      b.addDirichlet(0, 0);
    }
    return b;
  });

  /// Solve
  sshtSolver->solve();

  /// Print files
  petscVectopvtu(octDA, treePart,sshtSolver->getCurrentSolution(), "ssht", "ssht", varname, domainExtents, false, false, ndof);

  /// L2 error
  const auto analytic_sol = [](const TALYFEMLIB::ZEROPTV &pos, const DENDRITE_UINT dof, const DENDRITE_REAL time) {
#if (DIM == 2)
    return sin(M_PI * pos.x()) * sin(M_PI * pos.y());
#elif (DIM == 3)
    return sin(M_PI * pos.x()) * sin(M_PI * pos.y()) * sin(M_PI * pos.z());
#endif
  };
  VecInfo v(sshtSolver->getCurrentSolution(), 1, 0);
  Analytic sshtAnalytic(octDA,treePart, v, analytic_sol, physDomain);
  sshtAnalytic.getL2error();


  delete sshtEq;
  delete sshtSolver;
  dendrite_finalize(octDA);

}