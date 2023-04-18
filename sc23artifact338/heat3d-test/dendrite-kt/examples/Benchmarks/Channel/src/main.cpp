//
// Created by maksbh on 6/18/21.
//



#include <DendriteUtils.h>
#include <SubDA/SubDomain.h>
#include <feMatrix.h>
#include <IO/VTU.h>
#include <chrono>
#include <SDARefine.h>
int main(int argc, char *argv[]) {

  dendrite_init(argc, argv);
  if (argc < 4) {
    TALYFEMLIB::PrintStatus("Usage: ", argv[0], " level boundaryLevel eleOrder");
    return EXIT_FAILURE;
  }
  MPI_Barrier(MPI_COMM_WORLD);

  int level = std::atoi(argv[1]);
  int blevel = std::atoi(argv[2]);
  int eleOrder = std::atoi(argv[3]);

  DomainInfo fullDomain, physDomain;
  fullDomain.max.fill(8.0);
  fullDomain.min.fill(0.0);
  physDomain.max.fill(8.0);
  physDomain.min.fill(0.0);
  physDomain.max[1] = 1.0;
  physDomain.max[2] = 1.0;
  DomainExtents domain(fullDomain, physDomain);

  SubDomain subDomain(domain);
  std::function<ibm::Partition(const double *, double)> functionToRetain = [&](const double *physCoords,
                                                                               double physSize) {
    return (subDomain.functionToRetain(physCoords, physSize));
  };
  DistTREE dTree;
  DA *octDA = createSubDA(dTree, functionToRetain, level, eleOrder);
  subDomain.finalize(octDA,dTree.getTreePartFiltered(),domain);
  while (true) {
    Refine refine(octDA, dTree.getTreePartFiltered(), domain, blevel,level);
    DA *newDA = refine.getRefineSubDA(dTree);
    if (newDA == nullptr) {
      break;
    }
    std::swap(newDA, octDA);
    delete newDA;

    subDomain.finalize(octDA,dTree.getTreePartFiltered(),domain);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  TALYFEMLIB::PrintInfo("Started remeshing");
  auto start  = std::chrono::high_resolution_clock::now();
  for(int i = 0; i < 5; i++){
    Refine refine(octDA, dTree.getTreePartFiltered(), domain, blevel,level);
    DA *newDA = refine.getForceRefineSubDA(dTree);
    std::swap(octDA,newDA);
    delete newDA;
    TALYFEMLIB::PrintStatus("Global Node Sz = ", octDA->getGlobalNodeSz());
  }
  MPI_Barrier(MPI_COMM_WORLD);
  auto end  = std::chrono::high_resolution_clock::now();
  auto totalTime =  (static_cast<std::chrono::duration<DENDRITE_REAL>>(end - start)).count(); ;
  TALYFEMLIB::PrintStatus("TotalDA Creation Time = ", totalTime);

  MPI_Barrier(MPI_COMM_WORLD);
  dendrite_finalize(octDA);

}