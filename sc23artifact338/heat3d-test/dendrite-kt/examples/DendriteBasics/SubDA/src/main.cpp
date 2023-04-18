#include <iostream>
#include <distTree.h>
#include <oda.h>
#include <point.h>
#include <sfcTreeLoop_matvec_io.h>

#include <octUtils.h>
#include <PETSc/IO/petscVTU.h>
#include <DendriteUtils.h>
#include <SDARefine.h>
#include <PETSc/Solver/LinearSolver.h>
#include <PETSc/PetscUtils.h>
#include "SSHTEquation.h"
#include <SubDA/SubDomain.h>
#include <Boundary/SubDomainBoundary.h>

using namespace PETSc;
void checkHangingNodes(ot::DA<DIM>* octDA, const std::vector<TREENODE> & treePart, const DomainExtents & domainExtents){
  using ot::RankI;
  const size_t sz = octDA->getTotalNodalSz();
  auto partFront = octDA->getTreePartFront();
  auto partBack = octDA->getTreePartBack();
  const auto tnCoords = octDA->getTNCoords();
  double counter = 0;
  const std::vector<RankI> & ghostedGlobalNodeId = octDA->getNodeLocalToGlobalMap();
  std::ofstream fout("Hanging" + std::to_string(octDA->getRankActive()) +".txt");
  std::vector<double> hangingElements(octDA->getLocalElementSz(),0.0);

  ot::MatvecBaseIn<DIM,RankI,false> loop(sz,1,octDA->getElementOrder(), false,0,tnCoords, &(*ghostedGlobalNodeId.cbegin()),&(*treePart.cbegin()),treePart.size(),*partFront,*partBack);
  while(!loop.isFinished()){
    if (loop.isPre() && loop.subtreeInfo().isLeaf()) {
      if(loop.subtreeInfo().isElementBoundary()){
        if(true){
          const RankI * nodeIDs = loop.subtreeInfo().readNodeValsIn();
          const double *coords = loop.subtreeInfo().getNodeCoords();
          for(int i = 0; i < 8; i++){
            if(nodeIDs[i] == 0){
              hangingElements[counter] = 1.0;
              for(int j = 0; j < 8; j++) {
                fout << coords[DIM*j+0]*domainExtents.fullDADomain.max[0] << " " << coords[DIM*j+1]*domainExtents.fullDADomain.max[1] << " "<< coords[DIM*j+2]*domainExtents.fullDADomain.max[2] << "\n";
              }
              fout << "\n";
              std::cout << "Something wrong happened " << " "<< octDA->getRankActive() << "\n";
            }
          }


        }
      }
      counter++;
      loop.next();
    }
    else{
      loop.step();
    }
  }
  std::cout << counter << " "<< octDA->getLocalElementSz() << "\n";
  const char * varname[]{"hanging"};
  IO::writeVecTopVtu(octDA,treePart,hangingElements.data(),"Hanging","hanging",varname,domainExtents,true,false,1);
  fout.close();
}



int main(int argc, char * argv[]){
  dendrite_init(argc,argv);
  DENDRITE_UINT eleOrder = 2;
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
//  physDomain.max[0] = 1.0;
//  physDomain.max[1] = 1.0;
//  physDomain.max[0] = 300/64. ;
//  physDomain.max[1] = 150;
#endif
#if(DIM == 2)
  DomainInfo cubeDomain;
  cubeDomain.min.fill(0.0);
  cubeDomain.max.fill(2.0);
//  cubeDomain.max[0] = 1;
//  cubeDomain.max[1] = 1.;

  DomainInfo physDomain;
  physDomain.min.fill(0.0);
//  physDomain.max.fill(0.8);
  physDomain.max[0] = 2.0;
  physDomain.max[1] = 2.0;
#endif

  if(argc < 3){
    if(not(TALYFEMLIB::GetMPIRank())) {
      std::cout << "Usage: BaseRefinementLevel BlockRefinementLevel SphereRefinementLevel\n";
    }
    exit(EXIT_FAILURE);
  }

  DENDRITE_UINT baseLevel = std::atoi(argv[1]);
  DENDRITE_UINT blockLevel = std::atoi(argv[2]);
  DENDRITE_UINT objectLevel = std::atoi(argv[3]);

  OctToPhysical octToPhysical(cubeDomain);
  const Point<DIM> & scalingFactor = octToPhysical.getScalingFactor();

  DomainExtents domainExtents(cubeDomain, physDomain);
  SubDomain subDA(domainExtents);

#if(DIM == 3)
  double centerCoords[DIM] {physDomain.max[0]/(2),physDomain.max[1]/(2),physDomain.max[2]/(2)};
//  double _centerCoords[DIM];
//  static constexpr double RADIUS = 0.45;
//  static double centerSphere[3][3]{{0.5, 0.5, 8.0},
//                                   {0.5, 0.0, 8.0 + RADIUS * 2 + 0.01},
//                                   {0.5, 1.0, 8.0 + RADIUS * 2 + 0.01}};
//  _centerCoords[0] = centerSphere[0][0] ;
//  _centerCoords[1] = centerSphere[0][1] ;
//  _centerCoords[2] = centerSphere[0][2] ;
  subDA.addObject(VOXEL::Sphere(centerCoords,0.45));

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

  for(count = 0; count < 1; count++){
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
  auto sshtEq = new TalyEquation<SSHTEquation, SSHTNodeData>(octDA,treePart, domainExtents);
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
  petscVectopvtu(octDA, treePart,sshtSolver->getCurrentSolution(), "Solution","ssht", varname, domainExtents, false, false, ndof);
  dendrite_finalize(octDA);
#endif
#if(DIM == 2)
  std::function<ibm::Partition(const double *, double )> functionToRetain = [&](const double *physCoords, double physSize) {
    return (subDA.functionToRetain(physCoords, physSize));
  };

  DistTREE dtree;

  DA * octDA = createSubDA(dtree,functionToRetain,baseLevel,eleOrder);
  subDA.finalize(octDA,dtree.getTreePartFiltered(),domainExtents);
  IO::writeBoundaryElements(octDA,dtree.getTreePartFiltered(),"boundary","boundary",domainExtents);
  dendrite_finalize(octDA);
#endif

}
