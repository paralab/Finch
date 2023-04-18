//
// Created by maksbh on 5/29/21.
//
#include <DataTypes.h>
#include <DendriteUtils.h>
#include <Geometry/STL.h>
#include <SubDA/Voxel.h>
#include <SubDA/SubDomain.h>
#include <IO/VTU.h>
#include "Boundary/SubDomainBoundary.h"
#include <Refine.h>

void printStatistics(DA * octDA , const int refineLevel) {
  DENDRITE_UINT localElemSz = octDA->getLocalElementSz();
  DENDRITE_UINT globalElemSz;
  MPI_Reduce(&localElemSz, &globalElemSz,1, MPI_UINT32_T, MPI_SUM, 0, MPI_COMM_WORLD);
  TALYFEMLIB::PrintStatus(refineLevel," : Refinement level done");
  TALYFEMLIB::PrintStatus("Number of nodes: ", octDA->getGlobalNodeSz());
  TALYFEMLIB::PrintStatus("Number of elements: ", globalElemSz);
}
int main(int argc, char *argv[]) {

  dendrite_init(argc, argv);

  if (argc < 5) {
    if (TALYFEMLIB::GetMPIRank() == 0) {
      std::cout << "Usage " << argv[0] << "GeometryPath Level maxLevel IBM/CarvedOut[0,1]\n";
    }
    exit(EXIT_FAILURE);
  }

  const std::string stlFname = argv[1];
  const DENDRITE_UINT level = std::atoi(argv[2]);
  const DENDRITE_UINT maxLevel = std::atoi(argv[3]);
  RetainSide side = RetainSide::OUT;
  CaseType caseType = static_cast<CaseType>(std::atoi(argv[4]));

  DENDRITE_UINT eleOrder = 1;
  DENDRITE_UINT ndof = 1;
  m_uiMaxDepth = 25;

  DomainInfo cubeDomain;
  cubeDomain.min.fill(0.0);
  cubeDomain.max.fill(1.2);

  DomainExtents domainExtents(cubeDomain);
  OctToPhysical octToPhysical(domainExtents);

  std::array<DENDRITE_REAL, DIM> point;
  point.fill(0.0);

  GEOMETRY::STL stl(stlFname, GEOMETRY::InOutTest::RAY_TRACING);
  GEOMETRY::Geometry geometry(&stl, Point<DIM>(point), side);


  SubDomain subdomain(domainExtents);
  if (caseType == CaseType::CARVED) {
    subdomain.addObject(&geometry);
  }

  std::function<ibm::Partition(const double *, double)> functionToRetain = [&](const double *physCoords,
                                                                               double physSize) {
    return (subdomain.functionToRetain(physCoords, physSize));
  };

  DistTREE dtree;


  if (caseType == CARVED) {
    DA *octDA = createSubDA(dtree, functionToRetain, level, eleOrder);

    const auto &treeNode = dtree.getTreePartFiltered();
    subdomain.finalize(octDA, treeNode, domainExtents);
    for (int refineLevel = level; refineLevel <= maxLevel; refineLevel++) {
      while (true) {
        SubDomainBoundary subDomainBoundary(&subdomain, octDA, domainExtents);
        Refine sdaRefine(octDA, dtree.getTreePartFiltered(), domainExtents, &subDomainBoundary, refineLevel, caseType,
                         &stl,
                         false);
        DA *newDA = sdaRefine.getRefineSubDA(dtree);
        if (newDA == nullptr) {
          break;
        }
        std::swap(newDA, octDA);
        delete newDA;
        subdomain.finalize(octDA, dtree.getTreePartFiltered(), domainExtents);
      }
      printStatistics(octDA,refineLevel);
    }
  }
  if (caseType == IBM) {

    for (int refineLevel = level; refineLevel <= maxLevel; refineLevel++){
      DA *octDA = createSubDA(dtree, functionToRetain, level, eleOrder);
      while (true) {
        SubDomainBoundary subDomainBoundary(&subdomain, octDA, domainExtents);
        Refine sdaRefine(octDA, dtree.getTreePartFiltered(), domainExtents, &subDomainBoundary, refineLevel, caseType,&stl,false);
        DA *newDA = sdaRefine.getRefineSubDA(dtree);
        if (newDA == nullptr) {
          break;
        }
        std::swap(newDA, octDA);
        delete newDA;
        subdomain.finalize(octDA, dtree.getTreePartFiltered(), domainExtents);
      }
      SubDomainBoundary subDomainBoundary(&subdomain, octDA, domainExtents);
      Refine sdaRefine(octDA, dtree.getTreePartFiltered(), domainExtents, &subDomainBoundary, refineLevel, caseType,
                       &stl,
                       true);
      DA *newDA = sdaRefine.getRefineSubDA(dtree);
      if (newDA != NULL) {
        std::swap(newDA, octDA);
        delete newDA;
        subdomain.finalize(octDA, dtree.getTreePartFiltered(), domainExtents);
      }
      printStatistics(octDA,refineLevel);
      delete octDA;
    }
  }
}