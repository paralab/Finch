//
// Created by maksbh on 7/11/20.
//
#include <vector>
#include <DataTypes.h>
#include <point.h>
#include <DendriteUtils.h>
#include <IMGA/IMGA.h>
#include <SSHTNodeData.h>
#include <SSHTEquation.h>
#include <IMGA/Marker.h>
#include <TalyEquation.h>
#include <sfcTreeLoop_matvec_io.h>
#include <PETSc/PetscUtils.h>
#include <utils.h>
#include <PETSc/IO/petscVTU.h>
#include <CalcError.h>
#include <IMGALoop.h>
#include <IMGA/IMGAInterpolate.h>
#include <IMGA/MovingIMGA.h>
#include <IMGA/IO.h>

using namespace PETSc;

/**
 * @brief An example to test the background element of the processor.
 * In order for efficent implementation, Point should be used instead of ZEROPTV.
 */
int main(int argc, char *argv[]) {
  dendrite_init(argc, argv);

  GEOMETRY::STL stl(argv[1], GEOMETRY::InOutTest::RAY_TRACING);
  DomainInfo domainInfo;
  domainInfo.min.fill(0.0);
  domainInfo.max.fill(1.0);
  DomainExtents domainExtents(domainInfo);
  const DENDRITE_UINT eleOrder = 1;
  GEOMETRY::Geometry geometry(&stl, Point<DIM>(0,0,0));
  int level = 1;
  GeomRefinement geomRefine;
  geomRefine.maxSplitIteration = 0;
  geomRefine.octantLevel = level;
  geomRefine.ratioArea = 0.25;
  IMGA imga(domainExtents,&geometry,geomRefine);
  SubDomain subDomain(domainExtents);
  std::function<ibm::Partition(const double *, double )> functionToRetain = [&](const double *physCoords, double physSize) {
    return (subDomain.functionToRetain(physCoords, physSize));
  };
  DistTREE distTree;
  DA *octDA = createSubDA(distTree,functionToRetain,level,eleOrder);
  IO::writeSurfaceTopVTU<0>(octDA,&imga, nullptr,"surface","surface", nullptr);

#if(DIM == 2)
  SubDomain subDomain(domainExtents);
  double _center[DIM] {0.5,0.5};
  VOXEL::Circle circle(_center,0.5,RetainSide::IN);
  subDomain.addObject(circle);
  std::function<ibm::Partition(const double *, double )> functionToRetain = [&](const double *physCoords, double physSize) {
    return (subDomain.functionToRetain(physCoords, physSize));
  };
  DistTREE distTree;
  DA *octDA = createSubDA(distTree,functionToRetain,level,eleOrder);
  subDomain.finalize(octDA, distTree.getTreePartFiltered(), domainExtents);
  performRefinement(octDA,distTree,domainExtents,0,subDomain,0);
  IO::writeBoundaryElements(octDA,distTree.getTreePartFiltered(),"Boundary","boundary",domainExtents);
  const auto & treeNode = distTree.getTreePartFiltered();
  static constexpr DENDRITE_UINT numDivisions = 10000;
  static constexpr DENDRITE_REAL radius = 0.5;
  static constexpr DENDRITE_REAL dtheta = 2*M_PI/(numDivisions - 1);
  std::vector<std::array<DENDRITE_REAL,DIM>> points(numDivisions);
  for(int i = 0; i < numDivisions; i++){
    DENDRITE_REAL cosTheta = cos(i*dtheta);
    DENDRITE_REAL sinTheta = sin(i*dtheta);
    points[i][0] = radius*cosTheta;
    points[i][1] = radius*sinTheta;
  }
  GEOMETRY::MSH msh(points, points.size(),GEOMETRY::InOutTest2D::CIRCLE);
  static const DENDRITE_REAL center[DIM]{0.0,0.0};
  msh.initCircle(0.5,center);
  msh.correctNormals();
  Point<DIM> point(std::array<DENDRITE_REAL, DIM>{0.5, 0.5});
  GEOMETRY::Geometry geometry(&msh, point,RetainSide::IN);
  GeomRefinement geomRefine;
  geomRefine.maxSplitIteration = 0;
  geomRefine.octantLevel = level;
  geomRefine.ratioArea = 0.25;
  IMGA imga(domainExtents,&geometry,geomRefine,ibmMethod);
  imga.initIMGAComputation(octDA,treeNode);
//  printGaussPointsToFile("GaussPoints.txt",imga.getSurfaceGaussPoints());
  Marker marker(octDA,treeNode,domainExtents,&imga,MarkerType::GAUSS_POINT);
  marker.printMarker();
#endif

//  SubDomainBoundary subDomainBoundary(&subDomain,octDA,domainExtents);
//
//  auto sshtEq = new TalyEquation<SSHTEquation, SSHTNodeData>(octDA, treeNode, domainExtents,1, nullptr,(ibmMethod == NITSCHE)? false:true,
//                                                             (ibmMethod == NITSCHE)? nullptr:&subDomainBoundary,ibmMethod);
//  sshtEq->assignIBMConstructs(&imga,marker.getMarkers().data());
//  LinearSolver *sshtSolver = setLinearSolver(sshtEq, octDA, 1, false);
//  sshtSolver->solve();
//  static const char * varname[]{"u"};
//  petscVectopvtu(octDA,distTree.getTreePartFiltered(),sshtSolver->getCurrentSolution(),"Solution",varname,domainExtents,false,false,1);
//  const VecInfo vecInfo(sshtSolver->getCurrentSolution(), 1, 0, PLACEHOLDER_NONE);
//  CalcError calcError(octDA, distTree.getTreePartFiltered(), vecInfo, domainExtents, &subDomain);
//  double error[2];
//  calcError.getL2error(error);
//  const auto &elemError = calcError.getElementalError();
//  IO::writeVecTopVtu(octDA, distTree.getTreePartFiltered(), elemError.data(), "Error", "ElemError", varname,
//                     domainExtents, true);
//  IMGALoop loop(octDA,&imga,treeNode,vecInfo,domainExtents);
//  DENDRITE_REAL boundaryError[2];
//  loop.computeBoundaryError(boundaryError);
//  TALYFEMLIB::PrintStatus("Boundary Error =",boundaryError[0], " " , boundaryError[1]);
//  TALYFEMLIB::PrintStatus(error[0], " " , error[1]);
//  if (TALYFEMLIB::GetMPIRank() == 0) {
//    std::ofstream fout("Error.txt", std::ios::app);
//    fout << level << " " << error[0] << " " << error[1] << " " << boundaryError[0] << " " << boundaryError[1] << "\n";
//    fout.close();
//  }
//  MPI_Barrier(MPI_COMM_WORLD);
//
//
//
//  delete sshtEq;
//  delete sshtSolver;
dendrite_finalize(octDA);

}