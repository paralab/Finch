//
// Created by maksbh on 9/18/19.
//

#include <iostream>
#include <DendriteUtils.h>
#include <point.h>
#include <TalyEquation.h>
#include <TSHTNodeData.h>
#include <TSHTEquation.h>
#include <PETSc/Solver/LinearSolver.h>
#include <PETSc/PetscUtils.h>
#include <IO/VTU.h>
#include <PETSc/IO/petscVTU.h>
#include <PostProcessing/postProcessing.h>
#include <TSHTAnalytic.h>


using namespace PETSc;
int main(int argc, char * argv[]){

  dendrite_init(argc,argv);
  int rank = TALYFEMLIB::GetMPIRank();
  if (argc < 2) {
    if (not(rank))
      std::cout << "Usage: " << argv[0]
                << " eleOrder level"
                << std::endl;
    return 0;
  }
  TALYFEMLIB::PrintInfo("Total number of processor = ",TALYFEMLIB::GetMPISize());
  TALYFEMLIB::PrintInfo("size of DendroInt ", sizeof(DendroIntL));
  TALYFEMLIB::PrintInfo("size of PetscInt ", sizeof(PetscInt));

  const DENDRITE_UINT eleOrder = static_cast<DENDRITE_UINT>(std::atoi(argv[1]));
  const DENDRITE_UINT level = static_cast<DENDRITE_UINT>(std::atoi(argv[2]));
  TALYFEMLIB::PrintStatus("eleOrder ",eleOrder);
  TALYFEMLIB::PrintStatus("Level ",level);

  std::chrono::high_resolution_clock::time_point tDA_start = std::chrono::high_resolution_clock::now();
  ot::DA<DIM> * octDA = createRegularDA(level,eleOrder);
  std::chrono::high_resolution_clock::time_point tDA_end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> time_DA = std::chrono::duration_cast<std::chrono::duration<double>>(tDA_end - tDA_start);
  TALYFEMLIB::PrintStatus("DA-construction time = ",time_DA.count());
  int ndof = 1;
  bool mfree = true;
  std::array<double,DIM> Min;
  std::array<double,DIM> Max;

  Min.fill(0.0);
  Max.fill(1.0);

  Point<DIM> domainMin(Min);
  Point<DIM> domainMax(Max);



  auto talyEq = new TalyEquation<TSHTEquation,TSHTNodeData>(octDA,domainMin,domainMax);
  TALYFEMLIB::PrintStatus("Created DA");
  LinearSolver *solver = setLinearSolver(talyEq, octDA, ndof, mfree);

  solver->setBoundaryCondition([&](const TALYFEMLIB::ZEROPTV pos, unsigned int nodeID) -> Boundary {
    Boundary b;
    static constexpr double eps = 1e-14;
    double x = pos.x();
    double y = pos.y();
    double z = pos.z();


    bool on_wall = (fabs(x - domainMin.x(0)) < eps) ||
        (fabs(y - domainMin.x(1)) < eps) ||
        (fabs(z - domainMin.x(2)) < eps) ||
        (fabs(x - domainMax.x(0)) < eps) ||
        (fabs(y - domainMax.x(1)) < eps) ||
        (fabs(z - domainMax.x(2)) < eps);

#if (DIM == 4)
    double t = pos.t();
    on_wall = on_wall ||  (fabs(t) < eps);
#endif
    double val = exp(-t)*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z);
    if (on_wall) {
      b.addDirichlet(0, val);
    }
    return b;

  });
  MPI_Barrier(MPI_COMM_WORLD);
  std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
  solver->solve();
  MPI_Barrier(MPI_COMM_WORLD);
  std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
  TALYFEMLIB::PrintStatus("Total time  = ", time_span.count() );
//  petscVectopvtu(octDA,solver->getCurrentSolution(),"ssht", false, false,ndof);
  TSHTAnalytic sshtAnalytic(octDA,solver->getCurrentSolution(),ndof);
  TALYFEMLIB::PrintStatus("L2 error = ",sshtAnalytic.getL2Error());
  delete talyEq;
  delete solver;
  dendrite_finalize(octDA);


}