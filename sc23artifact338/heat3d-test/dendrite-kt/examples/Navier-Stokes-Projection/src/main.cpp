//
// Created by maksbh on 5/21/20.
//


#include <iostream>
#include <DendriteUtils.h>
#include <point.h>
#include <TalyEquation.h>
#include <NSNodeData.h>
#include <MomentumEquation.h>
#include <PETSc/Solver/LinearSolver.h>
#include <PETSc/PetscUtils.h>
#include <IO/VTU.h>
#include <Traversal/Analytic.h>
#include <PETSc/IO/petscVTU.h>
//#include <PostProcessing/postProcessing.h>
#include <NSUtils.h>
#include <NSBoundaryConditions.h>
#include <PressurePoissonEquation.h>
#include <VelocityUpdateEquation.h>

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

  NSInputData inputData;

  if (!(inputData.ReadFromFile())) {
    if (!rank) {
      throw TALYFEMLIB::TALYException() << "Can't read the config file \n";
    }
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  inputData.solverOptionsMomentum.apply_to_petsc_options("-momentum_");
  inputData.solverOptionsPressurePoisson.apply_to_petsc_options("-pp_");
  inputData.solverOptionsVelocityUpdate.apply_to_petsc_options("-vupdate_");
  NSParams nsParams(&inputData);

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

  DomainExtents domainExtents(inputData.physDomain);
  std::vector<TREENODE> treePart;
  ot::DA<DIM> *octDA = createRegularDA(treePart, level, eleOrder);

  /// perform initial refinement
  performRefinement(octDA,domainExtents, treePart,inputData);

  /// Time construct
  std::vector<DENDRITE_REAL> dt(1, 0.1);
  std::vector<DENDRITE_REAL> totalT(1, 1.0);
  TimeInfo ti(0.0, dt, totalT);

  /// Boundary condition
  NSBoundaryConditions nsBC(&inputData,&ti);

  /// Equation and solver setup
  auto momentumEq =
      new TalyEquation<MometumEquation, NSNodeData>(octDA,treePart, domainExtents, DIM, &ti,nullptr, &inputData,&nsParams);
  NonlinearSolver *momentumSolver = setNonLinearSolver(momentumEq, octDA, DIM, mfree);
  SNES mometum_snes = momentumSolver->snes();
  SNESSetOptionsPrefix(mometum_snes, "momentum_");
  SNESSetFromOptions(mometum_snes);

  auto pressurePoissonEq =
      new TalyEquation<PressurePoissonEquation, NSNodeData>(octDA, treePart, domainExtents, 1,&ti, nullptr, &inputData,&nsParams);
  LinearSolver * pressurePoissonSolver = setLinearSolver(pressurePoissonEq, octDA, 1, mfree);
  pressurePoissonEq->setTime(&ti);
  KSP pressurePoisson_ksp = pressurePoissonSolver->ksp();
  KSPSetOptionsPrefix(pressurePoisson_ksp, "pp_");
  KSPSetFromOptions(pressurePoisson_ksp);

  auto velocityUpdateEq =
      new TalyEquation<VelocityUpdateEquation, NSNodeData>(octDA,treePart,  domainExtents, DIM, &ti, nullptr, &inputData,&nsParams);
  LinearSolver * velocityUpdateSolver = setLinearSolver(velocityUpdateEq, octDA, DIM, mfree);
  velocityUpdateEq->setTime(&ti);
  KSP velocityUpdate_ksp = velocityUpdateSolver->ksp();
  KSPSetOptionsPrefix(velocityUpdate_ksp, "vupdate_");
  KSPSetFromOptions(velocityUpdate_ksp);

  /// Analytical solution to BC for MMS
  if (inputData.ifMMS) {
    nsBC.setAnalyticalFunction(AnalyticalSolution);
  }

  /// Boundary condition setup
  momentumSolver->setDirichletBoundaryCondition([&](const TALYFEMLIB::ZEROPTV &pos, int nodeID) -> Boundary {
    Boundary b;
    nsBC.getMomentumBoundaryCondition(b, pos);
    return b;
  });
  pressurePoissonSolver->setDirichletBoundaryCondition([&](const TALYFEMLIB::ZEROPTV &pos, int nodeID) -> Boundary {
    Boundary b;
    nsBC.getPressureBoundaryCondition(b, pos);
    return b;
  });
  velocityUpdateSolver->setDirichletBoundaryCondition([&](const TALYFEMLIB::ZEROPTV &pos, int nodeID) -> Boundary {
    Boundary b;
    return b;
  });



  //// Vectors
  Vec prev1VelocitySolution,prev2VelocitySolution;
  Vec prev1PressureSolution,prev2PressureSolution;

  octDA->petscCreateVector(prev1VelocitySolution, false, false, DIM);
  octDA->petscCreateVector(prev2VelocitySolution, false, false, DIM);
  octDA->petscCreateVector(prev1PressureSolution, false, false, 1);
  octDA->petscCreateVector(prev2PressureSolution, false, false, 1);
  VecSet(prev2VelocitySolution,0.0);
  VecSet(prev2PressureSolution,0.0);
  setInitialConditionVelocity(octDA, inputData,prev1VelocitySolution);
  setInitialConditionPressure(octDA, inputData,prev1PressureSolution);

  /// Sync vectors: For pressure poisson and velocity update matrix does not need vector information
  momentumEq->setVectors({
                         VecInfo(PLACEHOLDER_GUESS,    DIM, NSNodeData::VEL_X),
                         VecInfo(prev1VelocitySolution,DIM, NSNodeData::VEL_X_PRE1   ),
                         VecInfo(prev1PressureSolution,1,   NSNodeData::PRESSURE_PRE1),
                         VecInfo(prev2VelocitySolution,DIM, NSNodeData::VEL_X_PRE2   ),
                         VecInfo(prev2PressureSolution,1,   NSNodeData::PRESSURE_PRE2)
  },SYNC_TYPE::ALL);
  pressurePoissonEq->setVectors({
                        VecInfo(momentumSolver->getCurrentSolution(),DIM, NSNodeData::VEL_X),
                        VecInfo(prev1VelocitySolution,               DIM, NSNodeData::VEL_X_PRE1),
                        VecInfo(prev1PressureSolution,                1 , NSNodeData::PRESSURE_PRE1 ),
                        VecInfo(prev2VelocitySolution,               DIM, NSNodeData::VEL_X_PRE2    ),
                        VecInfo(prev2PressureSolution,                1 , NSNodeData::PRESSURE_PRE2 )
  },SYNC_TYPE::VECTOR_ONLY);
  velocityUpdateEq->setVectors({
                       VecInfo(momentumSolver->getCurrentSolution(),       DIM, NSNodeData::VEL_X        ),
                       VecInfo(pressurePoissonSolver->getCurrentSolution(), 1 , NSNodeData::PRESSURE     ),
                       VecInfo(prev1VelocitySolution,                      DIM, NSNodeData::VEL_X_PRE1   ),
                       VecInfo(prev1PressureSolution,                       1 , NSNodeData::PRESSURE_PRE1),
                       VecInfo(prev2VelocitySolution,                      DIM, NSNodeData::VEL_X_PRE2   ),
                       VecInfo(prev2PressureSolution,                       1 , NSNodeData::PRESSURE_PRE2)
  },SYNC_TYPE::VECTOR_ONLY);


  VecCopy(prev1VelocitySolution, momentumSolver->getCurrentSolution());
  VecCopy(prev1PressureSolution, pressurePoissonSolver->getCurrentSolution());

#if(DIM == 3)
  static const char *momentum_varname[]{"u", "v", "w"};
  static const char *pressure_varname[]{"p"};
#endif
#if(DIM == 2)
  static const char *momentum_varname[]{"u", "v"};
  static const char *pressure_varname[]{"p"};
#endif
  petscVectopvtu(octDA,treePart,  momentumSolver->getCurrentSolution(), "momentum_init", momentum_varname,domainExtents, false, false, DIM);
  petscVectopvtu(octDA,treePart, momentumSolver->getCurrentSolution(), "pressure_init", pressure_varname,domainExtents, false, false, 1);

  /// Setting first order for first step
  nsParams.updateNSOrder(1);

  while (ti.getCurrentTime() < ti.getEndTime()) {
    ti.increment();
    ti.print();

    /// Solve
    momentumSolver->solve();
    pressurePoissonSolver->solve();
    velocityUpdateSolver->solve();


    VecCopy(prev1VelocitySolution, prev2VelocitySolution);
    VecCopy(velocityUpdateSolver->getCurrentSolution(), prev1VelocitySolution);

    VecCopy(prev1PressureSolution, prev2PressureSolution);
    VecCopy(pressurePoissonSolver->getCurrentSolution(), prev1PressureSolution);

    nsParams.updateNSOrder(inputData.nsOrder);

    if(inputData.ifMMS) {
      {
        VecInfo v(prev1VelocitySolution, DIM, NSNodeData::VEL_X);
        Analytic NSAnalytic(octDA, treePart, v, AnalyticalSolution,domainExtents, ti.getCurrentTime());
        TALYFEMLIB::PrintStatus("Velocity error = ");
        NSAnalytic.getL2error();
      }
      {
        VecInfo v(prev1PressureSolution, 1 , NSNodeData::PRESSURE);
        Analytic NSAnalytic(octDA, treePart, v, AnalyticalSolution, domainExtents, ti.getCurrentTime());
        TALYFEMLIB::PrintStatus("Pressure error = ");
        NSAnalytic.getL2error();
      }

    }
  }

  petscVectopvtu(octDA, treePart,prev1VelocitySolution, "momentum", momentum_varname,domainExtents, false, false, DIM);
  petscVectopvtu(octDA, treePart, prev1PressureSolution, "pressure", pressure_varname,domainExtents, false, false, 1);

  delete momentumEq;
  delete momentumSolver;

  delete pressurePoissonEq;
  delete pressurePoissonSolver;

  delete velocityUpdateSolver;
  delete velocityUpdateEq;

  VecDestroy(&prev1VelocitySolution);
  VecDestroy(&prev2VelocitySolution);
  VecDestroy(&prev1PressureSolution);
  VecDestroy(&prev2PressureSolution);

  dendrite_finalize(octDA);

}