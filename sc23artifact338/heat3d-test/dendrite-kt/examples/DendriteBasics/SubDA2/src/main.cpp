//
// Created by maksbh on 11/8/20.
//

#include <DendriteUtils.h>
#include <SubDA/SubDomain.h>
#include <SDARefine.h>
#include <IO/VTU.h>
int main(int argc, char * argv[]){
    dendrite_init(argc,argv);
    DENDRITE_UINT eleOrder = 1;
    DENDRITE_UINT ndof = 1;
    m_uiMaxDepth = 25;

    MPI_Comm comm = MPI_COMM_WORLD;

    DomainInfo cubeDomain;
    cubeDomain.min.fill(0.00);
    cubeDomain.max.fill(26.64);

    DomainInfo physDomain;
    physDomain.min.fill(0.0);
    physDomain.max[0] = 15.0;
    physDomain.max[1] = 3.33;
    physDomain.max[2] = 26.64;


    if(argc < 3){
        if(not(TALYFEMLIB::GetMPIRank())) {
            std::cout << "Usage:  stlFname BaseLevel MaxLevel\n";
        }
        exit(EXIT_FAILURE);
    }


    std::string stlFname = argv[1];
    DENDRITE_UINT baseLevel = std::atoi(argv[2]);
    DENDRITE_UINT maxLevel = std::atoi(argv[3]);


    OctToPhysical octToPhysical(cubeDomain);
    const Point<DIM> & scalingFactor = octToPhysical.getScalingFactor();

    DomainExtents domainExtents(cubeDomain, physDomain);
    SubDomain subDA(domainExtents);


//#if(DIM == 2)
//    const DENDRITE_REAL center[DIM]{0.5,0.5};
//    VOXEL::Circle circle(center,0.4,RetainSide::IN);
//    subDA.addObject(circle);
//
//
//    const DENDRITE_REAL min[DIM]{0.4,0.4};
//    const DENDRITE_REAL max[DIM]{0.6,0.6};
//    VOXEL::Box box(min,max,RetainSide::OUT);
//    subDA.addObject(box);
//
//#endif
//#if(DIM == 3)
//    GEOMETRY::STL stl(stlFname);
//    std::array<DENDRITE_REAL,DIM> point; point.fill(0.0);
//    GEOMETRY::Geometry geometry(&stl,Point<DIM>(point),RetainSide::IN);
//
//
////    SubDomain subdomain(domainExtents);
//  subDA.addObject(&geometry);
////    const DENDRITE_REAL center[DIM]{0.5,0.5,0.5};
////    VOXEL::Sphere sphere(center,0.3,RetainSide::OUT);
////    subDA.addObject(sphere);
////    VOXEL::Sphere sphere1(center,0.1,RetainSide::OUT);
////    subDA.addObject(sphere1);
//#endif

    std::function<ibm::Partition(const double *, double )> functionToRetain = [&](const double *physCoords, double physSize) {
        return (subDA.functionToRetain(physCoords, physSize));
    };

    DistTREE dtree;

    DA * octDA = createSubDA(dtree,functionToRetain,baseLevel,eleOrder);
  subDA.finalize(octDA,dtree.getTreePartFiltered(),domainExtents);
  std::vector<ot::OCT_FLAGS::Refine> refineFlags(octDA->getLocalElementSz());
  std::fill(refineFlags.begin(), refineFlags.end(), ot::OCT_FLAGS::Refine::OCT_NO_CHANGE);
  DistTREE newDistTree, surrDistTree;
  DistTREE::distRemeshSubdomain(dtree, refineFlags, newDistTree, surrDistTree, 0.3);
  DA *newDA = new DA(newDistTree, MPI_COMM_WORLD, octDA->getElementOrder(), 100, 0.3);
  std::swap(newDistTree, dtree);
  std::swap(newDA, octDA);
  delete newDA;
  subDA.finalize(octDA,dtree.getTreePartFiltered(),domainExtents);

    IO::writeBoundaryElements(octDA,dtree.getTreePartFiltered(),"subDA","boundary",domainExtents);
    MPI_Barrier(MPI_COMM_WORLD);
    exit(EXIT_FAILURE);
    DENDRITE_UINT iter = 0;
    while (true){
        SDARefine sdaRefine(octDA,dtree.getTreePartFiltered(),domainExtents,maxLevel);
        DA* newDA = sdaRefine.getRefineSubDA(dtree);
        if(newDA == nullptr){
          break;
        }
        std::swap(newDA,octDA);
        delete newDA;
        subDA.finalize(octDA,dtree.getTreePartFiltered(),domainExtents);
        TALYFEMLIB::PrintStatus(iter++ , "iteration of Refinement done. Number of nodes:", octDA->getGlobalNodeSz());
    }

    IO::writeBoundaryElements(octDA,dtree.getTreePartFiltered(),"subDARefined","boundary",domainExtents);
    dendrite_finalize(octDA);

}