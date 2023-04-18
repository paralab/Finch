//
// Created by maksbh on 1/18/21.
//

#include <DataTypes.h>
#include <DendriteUtils.h>
#include <Geometry/STL.h>
#include <SubDA/Voxel.h>
#include <SubDA/SubDomain.h>
#include <IO/VTU.h>
#include <SDARefine.h>
#include "Boundary/SubDomainBoundary.h"

void printBoundaryNodesIndices(DA *octDA, const DENDRITE_UINT iter, const OctToPhysical &octToPhysical) {
  std::string foldername = "refineLevel" + std::to_string(iter);
  if (not(TALYFEMLIB::GetMPIRank())) {
    int ierr = mkdir(foldername.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    if (ierr != 0 && errno != EEXIST) {
      TALYFEMLIB::PrintError("Could not create folder for storing results (", strerror(errno), ").");
      return;
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  std::vector<std::size_t> bdy_index;
  octDA->getBoundaryNodeIndices(bdy_index);
  double coords[DIM];
  std::string fname = foldername + "/rank" + std::to_string(octDA->getRankActive()) + ".txt";
  std::ofstream fout(fname);
  // Looping over the outer boundaries to check if the dirichlet BC needs to be specified
  const std::vector<DendroIntL> &localToGlobalMap = octDA->getNodeLocalToGlobalMap();
  for (std::size_t bdy: bdy_index) {
    ot::treeNode2Physical(octDA->getTNCoords()[bdy + octDA->getLocalNodeBegin()], octDA->getElementOrder(), coords);
    TALYFEMLIB::ZEROPTV zeroptvCoords;
    octToPhysical.convertCoordsToPhys(zeroptvCoords, coords);
    if (not((zeroptvCoords.x() == 0) or (zeroptvCoords.x() == 1) or (zeroptvCoords.y() == 0) or (zeroptvCoords.y() == 1)
        or (zeroptvCoords.z() == 0) or (zeroptvCoords.z() == 1)))
      fout << zeroptvCoords[0] << " " << zeroptvCoords[1] << " " << zeroptvCoords[2] << "\n";
  }
  fout.close();
}
int main(int argc, char *argv[]) {
  dendrite_init(argc, argv);

  if (argc < 5) {
    if (TALYFEMLIB::GetMPIRank() == 0) {
      std::cout << "Usage " << argv[0] << "GeometryPath Level maxLevel RetainSide\n";
    }
    exit(EXIT_FAILURE);
  }

  const std::string stlFname = argv[1];
  const DENDRITE_UINT level = std::atoi(argv[2]);
  const DENDRITE_UINT maxLevel = std::atoi(argv[3]);
  RetainSide side = static_cast<RetainSide>(std::atoi(argv[4]));

  DENDRITE_UINT eleOrder = 1;
  DENDRITE_UINT ndof = 1;
  m_uiMaxDepth = 25;

  DomainInfo cubeDomain;
  cubeDomain.min.fill(0.0);
  cubeDomain.max.fill(1.0);

  DomainExtents domainExtents(cubeDomain);
  OctToPhysical octToPhysical(domainExtents);

  std::array<DENDRITE_REAL, DIM> point;
  point.fill(0.0);
#if (DIM == 3)
  GEOMETRY::STL stl(stlFname, GEOMETRY::InOutTest::RAY_TRACING);
  GEOMETRY::Geometry geometry(&stl, Point<DIM>(point), side);
#endif
#if (DIM == 2)
  GEOMETRY::MSH msh(stlFname, GEOMETRY::InOutTest2D::RAY_TRACING_2D);
  GEOMETRY::Geometry geometry(&msh, Point<DIM>(point), side);
#endif
  SubDomain subdomain(domainExtents);
  subdomain.addObject(&geometry);

  std::function<ibm::Partition(const double *, double)> functionToRetain = [&](const double *physCoords, double physSize) {
    return (subdomain.functionToRetain(physCoords, physSize));
  };

  DistTREE dtree;

  DA *octDA = createSubDA(dtree, functionToRetain, level, eleOrder);

  const auto &treeNode = dtree.getTreePartFiltered();
  subdomain.finalize(octDA, treeNode, domainExtents);
  printBoundaryNodesIndices(octDA, 0, octToPhysical);

  for(int refineLevel = level; refineLevel <= maxLevel; refineLevel++){
  while (true) {
    SubDomainBoundary subDomainBoundary(&subdomain, octDA, domainExtents);
    SDARefine sdaRefine(octDA, dtree.getTreePartFiltered(), domainExtents, &subDomainBoundary, refineLevel);
    DA *newDA = sdaRefine.getRefineSubDA(dtree);
    if (newDA == nullptr) {
      break;
    }
    std::swap(newDA, octDA);
    delete newDA;
    subdomain.finalize(octDA, dtree.getTreePartFiltered(), domainExtents);
  }
    printBoundaryNodesIndices(octDA, refineLevel, octToPhysical);
    DENDRITE_UINT  localElemSz = octDA->getLocalElementSz();
    DENDRITE_UINT  globalElemSz;
    MPI_Reduce(&localElemSz,&globalElemSz,1,MPI_UINT32_T,MPI_SUM,0,MPI_COMM_WORLD);
    TALYFEMLIB::PrintStatus(refineLevel, " : Refinement level done");
    TALYFEMLIB::PrintStatus("Number of nodes:", octDA->getGlobalNodeSz());
    TALYFEMLIB::PrintStatus("Number of elements:", globalElemSz);
    IO::writeBoundaryElements(octDA, treeNode, ("boundary"+std::to_string(refineLevel)).c_str(), "boundary", domainExtents);
  }



  MPI_Barrier(MPI_COMM_WORLD);

  dendrite_finalize(octDA);

}