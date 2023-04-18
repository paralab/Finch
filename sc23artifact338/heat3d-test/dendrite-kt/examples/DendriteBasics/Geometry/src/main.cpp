//
// Created by maksbh on 8/12/20.
//

#include <DataTypes.h>
#include <DendriteUtils.h>
#include <Geometry/STL.h>
#include <SubDA/Voxel.h>
#include <SubDA/SubDomain.h>
#include <IO/VTU.h>
#include <SDARefine.h>
int main(int argc, char * argv[]) {
  dendrite_init(argc, argv);
  if(argc < 2){
    if(not (TALYFEMLIB::GetMPIRank())){
      std::cout << "Usage : " << argv[0] << " STLPath maxRefineLevel" << "\n";
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
  }
  DENDRITE_UINT eleOrder = 1;
  DENDRITE_UINT ndof = 1;
  m_uiMaxDepth = 25;

  DomainInfo cubeDomain;
  cubeDomain.min.fill(0.0);
  cubeDomain.max.fill(1.2);
  DENDRITE_UINT  initialLevel = 3;

//
//
  DomainExtents domainExtents(cubeDomain);
//
//  std::cout << "Reading\n";
  std::array<DENDRITE_REAL,DIM> point; point.fill(0.0);
#if (DIM == 3)
  GEOMETRY::STL stl(argv[1], GEOMETRY::InOutTest::RAY_TRACING);
  GEOMETRY::Geometry geometry(&stl,Point<DIM>(point),RetainSide::OUT);
  DENDRITE_UINT  maxLevel = std::atoi(argv[2]);
#endif
#if (DIM == 2)
  GEOMETRY::MSH msh(argv[1], GEOMETRY::InOutTest2D::RAY_TRACING_2D);
  GEOMETRY::Geometry geometry(&msh,Point<DIM>(point),RetainSide::IN);
#endif
////  point.fill(2.0);
////  geometry.addTranslations(Point<DIM>(point));
//
//
  SubDomain subdomain(domainExtents);
  subdomain.addObject(&geometry);
//#if (DIM == 3)
////  std::cout << subdomain.getGeometry()[0]->getSTL()->getTriangles().data() << " "<< stl.getTriangles().data() << "\n";
//#endif
//#if (DIM == 2)
//  std::cout << subdomain.getGeometry()[0]->getMSH()->getLines().data() << " "<< msh.getLines().data() << "\n";
//#endif
  std::function<ibm::Partition(const double *, double )> functionToRetain = [&](const double *physCoords, double physSize) {
    return (subdomain.functionToRetain(physCoords, physSize));
  };
//
  DistTREE dtree;
//
  DA * octDA = createSubDA(dtree,functionToRetain,initialLevel,eleOrder);
  const auto & treeNode = dtree.getTreePartFiltered();
  subdomain.finalize(octDA, treeNode,domainExtents);
  for(int i = initialLevel; i < maxLevel; i++){
    while (true) {
      SubDomainBoundary subDomainBoundary(&subdomain, octDA, domainExtents);
      SDARefine sdaRefine(octDA, dtree.getTreePartFiltered(), domainExtents, &subDomainBoundary, i);
      DA * newDA = sdaRefine.getRefineSubDA(dtree);
      if(newDA == nullptr){
        break;
      }
      std::swap(octDA,newDA);
      delete newDA;
    }

    IO::writeBoundaryElements(octDA,treeNode,("boundary_lev"+std::to_string(i)).c_str(),"boundary",domainExtents);
    MPI_Barrier(MPI_COMM_WORLD);
  }




  dendrite_finalize(octDA);


}