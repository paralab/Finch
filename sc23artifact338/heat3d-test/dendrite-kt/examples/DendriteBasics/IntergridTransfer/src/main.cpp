//
// Created by maksbh on 11/9/20.
//

#include <DendriteUtils.h>
#include <DARefine.h>
#include <PETSc/IO/petscVTU.h>
#include <SubDA/SubDomain.h>
#include <DARefineValues.h>
using namespace PETSc;
static const char *varname[]{"u","v"};


int main(int argc, char * argv[]){
    dendrite_init(argc,argv);

    if(argc < 3){
        std::cout << "Usage: " << argv[0] << " "<< "level eleOrder\n";
        exit(EXIT_FAILURE);
    }

    /// Domain Extents
    DomainInfo cubeDomain, physDomain;
    cubeDomain.min.fill(0.0);
    cubeDomain.max.fill(4.0);
    physDomain.min.fill(0.0);
    physDomain.max.fill(4.0);
    physDomain.max[0] = 4.0;
    DomainExtents domainExtents(cubeDomain,physDomain);

    const DENDRITE_UINT level = static_cast<DENDRITE_UINT>(std::atoi(argv[1]));
    const DENDRITE_UINT eleOrder = static_cast<DENDRITE_UINT>(std::atoi(argv[2]));

    /// SubDomain Objects
    DistTREE distTree;
    SubDomain subDomain(domainExtents);
    std::function<ibm::Partition(const double *, double )> functionToRetain = [&](const double *physCoords, double physSize) {
        return (subDomain.functionToRetain(physCoords, physSize));
    };

    /// Create DA
    DA * octDA = createSubDA(distTree,functionToRetain,level,eleOrder,0.3);
    std::cout << octDA->getLocalNodalSz() << "\n";

    IO::writeBoundaryElements(octDA,distTree.getTreePartFiltered(),"Boundary","bnd",domainExtents);
    Vec initialVector;
    octDA->petscCreateVector(initialVector,false,false,1);
    OctToPhysical octToPhysical(domainExtents);
    int counter = 0;
    std::function<void(const double *, double *)> initial_condition = [&](const double *x, double *var) {
      TALYFEMLIB::ZEROPTV coords;
      std::memcpy(coords.data(),x, sizeof(double )*DIM);
      octToPhysical.convertCoordsToPhys(coords.data());
      var[0] = coords[0] + coords[1];
    };
    octDA->petscSetVectorByFunction(initialVector,initial_condition,false,false,1);
    static const char * varname[]{"phi"};
    petscVectopvtu(octDA,distTree.getTreePartFiltered(),initialVector,"initial","initial",varname,domainExtents,false,false,1);

//  while(true){
//    DARefine daRefine(octDA, distTree.getTreePartFiltered(), domainExtents,1);
//    DA *newDA = daRefine.getRefineSubDA(distTree);
//    if (newDA == nullptr) {
//      break;
//    }
//    daRefine.petscIntergridTransfer(newDA, distTree, initialVector, 1);
//    std::swap(octDA, newDA);
//    delete newDA;
//  }
//  petscVectopvtu(octDA, distTree.getTreePartFiltered(), initialVector, "newDA", "newDA", varname, domainExtents, false,
//                 false, 1);
  {
    DARefine daRefine(octDA, distTree.getTreePartFiltered(), domainExtents, 1, true);
    DA *newDA = daRefine.getRefineSubDA(distTree);
    daRefine.initPetscIntergridTransfer();
    daRefine.petscIntergridTransfer(newDA, distTree, initialVector, 1);
    daRefine.finializeIntergridTransfer();
    std::swap(octDA, newDA);
    delete newDA;
  }
  petscVectopvtu(octDA,distTree.getTreePartFiltered(),initialVector,"initial","final",varname,domainExtents,false,false,1);
////      /// Intergrid Transfer loop
//      int i = 0;
//      counter = 1;
//      while (true) {
//        i++;
//        octDA->petscSetVectorByFunction(initialVector,initial_condition,false,false,1);
//        TALYFEMLIB::PrintStatus("Active comm = ", octDA->getNpesActive());
//        DARefineValues daRefine(octDA, distTree.getTreePartFiltered(),initialVector, domainExtents, STAGE::REFINE);
//
//        DA *newDA = daRefine.getRefineSubDA(distTree);
//        if (newDA == nullptr) {
//          break;
//        }
//        daRefine.petscIntergridTransfer(newDA, distTree, initialVector, 1);
//        std::swap(octDA, newDA);
//        delete newDA;
//        petscVectopvtu(octDA, distTree.getTreePartFiltered(), initialVector, "newDA", "newDA", varname, domainExtents, false,
//                       false, 1);
//      }
//      TALYFEMLIB::PrintStatus("Coarsening stage");
//  DARefineValues daRefine(octDA, distTree.getTreePartFiltered(),initialVector, domainExtents, STAGE::COARSE);
//
//  DA *newDA = daRefine.getRefineSubDA(distTree);
//
//  daRefine.petscIntergridTransfer(newDA, distTree, initialVector, 1);
//  std::swap(octDA, newDA);
//  delete newDA;
//  petscVectopvtu(octDA, distTree.getTreePartFiltered(), initialVector, "newDA", "newDA", varname, domainExtents, false,
//                 false, 1);


//  for(int count = 1; count < 3; count++) {
//    int i = 0;
//    while (true) {
//      i++;
//      TALYFEMLIB::PrintStatus("Active comm = ", octDA->getNpesActive());
//      DARefine daRefine(octDA, distTree.getTreePartFiltered(), domainExtents, true,count);
////      daRefine.printefineFlags(("refineFlags" + std::to_string(i)).c_str(), "refine", domainExtents);
//
//      DA *newDA = daRefine.getRefineSubDA(distTree);
//      if (newDA == nullptr) {
//        break;
//      }
//      daRefine.petscIntergridTransfer(newDA, distTree, oldDAVec, ndof);
//      std::swap(octDA, newDA);
//      delete newDA;
//      petscVectopvtu(octDA, distTree.getTreePartFiltered(), oldDAVec, ("newDA"+std::to_string(count)).c_str(), "newDA", varname, domainExtents, false,
//                     false, ndof);
//    }
//      for (int i = 0; i < 1; i++) {
//        std::cout << "Coarsening start " << octDA->getLocalElementSz() << "\n";
//        DARefine daRefine(octDA, distTree.getTreePartFiltered(), domainExtents, false,count);
//        daRefine.printefineFlags("refineCoarsen", "coarsen", domainExtents);
//        DA *newDA = daRefine.getRefineSubDA(distTree);
//
//        daRefine.petscIntergridTransfer(newDA, distTree, oldDAVec, ndof);
//        std::swap(octDA, newDA);
//        delete newDA;
//        petscVectopvtu(octDA, distTree.getTreePartFiltered(), oldDAVec, ("newDACoarsened"+std::to_string(count)).c_str(), "newDA", varname,
//                       domainExtents,
//                       false,
//                       false, ndof);
//        std::cout << "Coarsening complete " << octDA->getLocalElementSz() << "\n";
//      }
//    }
//    /// Bunch of stuff to delete
    dendrite_finalize(octDA);
}