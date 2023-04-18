//
// Created by maksbh on 10/15/20.
//
#include <DendriteUtils.h>
#include <SubDA/SubDomain.h>
#include <Boundary/SubDomainBoundary.h>
#include "Loop.h"
#include <IO/VTU.h>
#include <PETSc/VecBounds.h>
int main(int argc, char * argv[]){
    dendrite_init(argc,argv);
    DENDRITE_UINT eleOrder = 1;
    DENDRITE_UINT ndof = 1;
    m_uiMaxDepth = 25;

    MPI_Comm comm = MPI_COMM_WORLD;
#if(DIM == 3)
    DomainInfo cubeDomain;
  cubeDomain.min.fill(0.0);
  cubeDomain.max.fill(8.0);
//  cubeDomain.max[1] = 1.0;
//  cubeDomain.max[0] = 1.0;

  DomainInfo physDomain;
  physDomain.min.fill(0.0);
  physDomain.max.fill(8.0);
  physDomain.max[0] = 1.0;
  physDomain.max[1] = 1.0;
//  physDomain.max[0] = 300/64. ;
//  physDomain.max[1] = 150;
#endif
#if(DIM == 2)
    DomainInfo cubeDomain;
  cubeDomain.min.fill(0.0);
  cubeDomain.max.fill(1.0);
//  cubeDomain.max[0] = 1;
//  cubeDomain.max[1] = 1.;

  DomainInfo physDomain;
  physDomain.min.fill(0.0);
//  physDomain.max.fill(0.8);
  physDomain.max[0] = 1.0;
  physDomain.max[1] = 0.5;
#endif

    if(argc < 2){
        if(not(TALYFEMLIB::GetMPIRank())) {
            std::cout << "Usage: Level\n";
        }
        exit(EXIT_FAILURE);
    }

    DENDRITE_UINT baseLevel = std::atoi(argv[1]);


    OctToPhysical octToPhysical(cubeDomain);
    const Point<DIM> & scalingFactor = octToPhysical.getScalingFactor();

    DomainExtents domainExtents(cubeDomain, physDomain);
    SubDomain subDA(domainExtents);

#if(DIM == 3)
    //  double centerCoords[DIM] {physDomain.max[0]/(2),physDomain.max[1]/(2),physDomain.max[2]/(2)};
//  double _centerCoords[DIM];
//  static constexpr double RADIUS = 0.45;
//  static double centerSphere[3][3]{{0.5, 0.5, 8.0},
//                                   {0.5, 0.0, 8.0 + RADIUS * 2 + 0.01},
//                                   {0.5, 1.0, 8.0 + RADIUS * 2 + 0.01}};
//  _centerCoords[0] = centerSphere[0][0] ;
//  _centerCoords[1] = centerSphere[0][1] ;
//  _centerCoords[2] = centerSphere[0][2] ;
//  subDA.addObject(VOXEL::Sphere(_centerCoords,0.45));
//
//  _centerCoords[0] = centerSphere[1][0] ;
//  _centerCoords[1] = centerSphere[1][1] ;
//  _centerCoords[2] = centerSphere[1][2] ;
//  subDA.addObject(VOXEL::Sphere(_centerCoords,0.45));
//
//  _centerCoords[0] = centerSphere[2][0] ;
//  _centerCoords[1] = centerSphere[2][1] ;
//  _centerCoords[2] = centerSphere[2][2] ;
//  subDA.addObject(VOXEL::Sphere(_centerCoords,0.45));

//  for(int l = 0; l < 2; l++){
//  for(int i = 0; i < 2; i++) {
//
//    _centerCoords[0] = centerCoords[0];
//    _centerCoords[1] = centerCoords[1] - 1.2 * i ;//+ (l%2)*0.6;
//    _centerCoords[2] = centerCoords[2] + l*1.2;
//    if (_centerCoords[1] - 0.5 > 0) {
//      subDA.addObject(VOXEL::Sphere(_centerCoords, 0.45));
//    }
//
//
//    _centerCoords[0] = centerCoords[0];
//    _centerCoords[1] = centerCoords[1] + 1.2 * i ;//+ (l%2)*0.6;
//    _centerCoords[2] = centerCoords[2] + l*1.2;
//    if (_centerCoords[1] + 0.5 < physDomain.max[1]) {
//      subDA.addObject(VOXEL::Sphere(_centerCoords, 0.45));
//    }
//  }
//
//
//  }

  std::function<ibm::Partition(const double *, double )> functionToRetain = [&](const double *physCoords, double physSize) {
    return (subDA.functionToRetain(physCoords, physSize));
  };

  DistTREE distTree;

  DA * octDA = createSubDA(distTree,functionToRetain,baseLevel,eleOrder);
  const auto & treePart = distTree.getTreePartFiltered();
  subDA.finalize(octDA,treePart,domainExtents);
  std::cout << "Num Local Element = " << octDA->getLocalElementSz() << "\n";
  IO::writeBoundaryElements(octDA,treePart,"subDA","subDA",domainExtents);
  MPI_Barrier(MPI_COMM_WORLD);


  DENDRITE_UINT count = 0;

  for(count = 0; count < 4; count++){
    SDARefine refine(octDA, treePart,domainExtents,subDA, blockLevel, objectLevel);
    DA * newDA = refine.getRefineSubDA(distTree);
    if(newDA == nullptr){
      break;
    }
    std::swap(octDA, newDA);
    delete newDA;

    TALYFEMLIB::PrintInfo(count, " Iteration complete"," ",treePart.size()," ", octDA->getLocalElementSz());
    IO::writeBoundaryElements(octDA,treePart,("RefineSubDA"+std::to_string(count)).c_str(),"subDA",domainExtents);

    count++;
      subDA.finalize(octDA,treePart,domainExtents);

  }

  IO::writeBoundaryElements(octDA, distTree.getTreePartFiltered(), "BoundaryRefined", "Boundary", domainExtents);
  checkHangingNodes(octDA,distTree.getTreePartFiltered(),domainExtents);
//  checkHangingNodes(octDA,treePart,domainExtents);
  /*auto sshtEq = new TalyEquation<SSHTEquation, SSHTNodeData>(octDA, domainExtents);
  LinearSolver *sshtSolver = setLinearSolver(sshtEq, octDA, ndof, false);

  /// Boundary condition
  sshtSolver->setDirichletBoundaryCondition([&](const TALYFEMLIB::ZEROPTV pos, unsigned int nodeID) -> Boundary {
    Boundary b;
    b.addDirichlet(0, 10);
    return b;
  });

  /// Solve
  sshtSolver->solve();

  const char*varname[]{"T"};
  /// Print files
  petscVectopvtu(octDA, sshtSolver->getCurrentSolution(), "Solution","ssht", varname, domainExtents, false, false, ndof);*/
  dendrite_finalize(octDA);
#endif
#if(DIM == 2)
    std::function<ibm::Partition(const double *, double )> functionToRetain = [&](const double *physCoords, double physSize) {
    return (subDA.functionToRetain(physCoords, physSize));
  };

  DistTREE dtree;

  DA * octDA = createSubDA(dtree,functionToRetain,baseLevel,eleOrder);
  subDA.finalize(octDA,dtree.getTreePartFiltered(),domainExtents);
  IO::writeBoundaryElements(octDA,dtree.getTreePartFiltered(),"subDA","boundary",domainExtents);
  SubDomainBoundary subDomainBoundary(&subDA,octDA,domainExtents);
    Vec localVec;
    octDA->petscCreateVector(localVec,false,false,1);
    VecInfo vec(localVec,1,0);
  Loop loop(octDA,dtree.getTreePartFiltered(),vec,domainExtents,&subDomainBoundary);
  loop.getTotalValue();

  Vec vecTest;
  octDA->petscCreateVector(vecTest,false,false,2);
  std::function<void(const double *, double *)> initial_condition = [](const double *x, double *var) {
    var[0] = 1;
    var[1] = 2;
  };
  octDA->petscSetVectorByFunction(vecTest,initial_condition,false,false,2);
  PetscScalar max[2],min[2];
  VecBounds<2>::getMaxAndMinimumValues(vecTest,octDA->getCommActive(),max,min);



//  vecBounds.getMaxAndMinimumValues(max,min);
  std::cout << max[0] << " " << max[1] << "\n";
  std::cout << min[0] << " " << min[1] << "\n";
  max[0] = 1; max[1] = 1;
  min[0] = 1; min[1] = 1;
  VecBounds<2>::setMaxAndMinimumValues(vecTest,octDA->getCommActive(),max,min);
  VecBounds<2>::getMaxAndMinimumValues(vecTest,octDA->getCommActive(),max,min);
  std::cout << max[0] << " " << max[1] << "\n";
  std::cout << min[0] << " " << min[1] << "\n";

  dendrite_finalize(octDA);

#endif

}