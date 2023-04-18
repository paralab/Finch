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
#include <PostProcessing/postProcessing.h>
#include <SSHTAnalytic.h>

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
  const DENDRITE_UINT eleOrder = static_cast<DENDRITE_UINT>(std::atoi(argv[1]));
  const DENDRITE_UINT level = static_cast<DENDRITE_UINT>(std::atoi(argv[2]));
  TALYFEMLIB::PrintStatus("eleOrder ",eleOrder);
  TALYFEMLIB::PrintStatus("Level ",level);
  TALYFEMLIB::PrintStatus("DIM ",DIM);
  ot::DA<DIM> subDA;


  std::array<DENDRITE_UINT,DIM> extend;

  std::vector<TREENODE> treeNode;
#if(DIM == 3)
  extend = {level-2,level-1,level-2};
#elif(DIM == 2)
  extend = {level-2,level-1};
#endif
  createSubDA(subDA,treeNode,extend,level,eleOrder);
  std::vector<ot::OCT_FLAGS::Refine> refineFlags(treeNode.size(),ot::OCT_FLAGS::Refine::OCT_NO_CHANGE);
  refineFlags[0] = ot::OCT_FLAGS::Refine::OCT_REFINE;
  std::vector<ot::TreeNode<DENDRITE_UINT, DIM>>  newTree;
  std::vector<ot::TreeNode<DENDRITE_UINT, DIM>>  surrTree;
  ot::SFC_Tree<DENDRITE_UINT , DIM>::distRemesh(treeNode, refineFlags, newTree, surrTree, 0.3, MPI_COMM_WORLD);
  DA * newDA = new ot::DA<DIM>(newTree,MPI_COMM_WORLD,eleOrder,100,0.3);
  std::cout << newTree.size() << "\n";
//  subDA = *newDA;
//  delete subDA;

  int ndof = 1;
  bool mfree = true;

  Point<DIM> domainMin(0.0,0.0,0.0);
  Point<DIM> domainMax(1.0,1.0,1.0);


  auto talyEq = new TalyEquation<SSHTEquation,SSHTNodeData>(newDA,domainMin,domainMax);
  LinearSolver *solver = setLinearSolver(talyEq, newDA, ndof, mfree);

  solver->setDirichletBoundaryCondition([&](const TALYFEMLIB::ZEROPTV pos, unsigned int nodeID) -> Boundary {
    Boundary b;
    b.addDirichlet(0, 0);
    return b;
  });

  solver->solve();
  const char *varname[]{"u"};
//  petscVectopvtu(newDA,solver->getCurrentSolution(),"ssht",varname, false, false,ndof);
  VecInfo vec(solver->getCurrentSolution(),ndof,0);

  TSHTAnalytic sshtAnalytic(newDA,vec);
  TALYFEMLIB::PrintStatus("L2 error = ",sshtAnalytic.getL2Error());
  MPI_Barrier(MPI_COMM_WORLD);
  delete talyEq;
  delete solver;
  dendrite_finalize();


}