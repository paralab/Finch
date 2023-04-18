
// Created by maksbh on 7/27/20.


#include <Checkpoint/Checkpointer.h>
#include <DendriteUtils.h>
#include <SubDA/SubDomain.h>
#include <IO/VTU.h>
#include <PETSc/IO/petscVTU.h>
#include <SDARefine.h>
static const char *varname[]{"u","v"};
int main(int argc, char *argv[]){

  dendrite_init(argc,argv);
  int eleOrder = 1;
  int level = 4 ;
  bool resumeFromCheckpoint = false;
  if(argc > 1){
    resumeFromCheckpoint = static_cast<bool>(std::atoi(argv[1]));
  }
  static constexpr int ndof = 2;

//  TALYFEMLIB::
  DomainInfo fullDomain,physicalDomain;
  fullDomain.max.fill(8.0);
  fullDomain.min.fill(0.0);
  physicalDomain.max.fill(8.0);
  physicalDomain.max[1] = 2.0;
  physicalDomain.min.fill(0.0);
  DomainExtents domain(fullDomain,physicalDomain);

  SubDomain subDA(domain,resumeFromCheckpoint);
  double center[3]{4.5,0.6,4.5};
  VOXEL::Sphere sphere(center,0.5,RetainSide::OUT);
  subDA.addObject(sphere);
  std::function<ibm::Partition(const double *, double )> functionToRetain = [&](const double *physCoords, double physSize) {
    return (subDA.functionToRetain(physCoords, physSize));
  };
  DA *octDA;
  DistTREE dtree;
  if(resumeFromCheckpoint){
    TALYFEMLIB::PrintStatus("Loading checkpoint");
    Checkpointer checkpointer(5, "Checkpoint");
    std::vector<VecInfo> vecs;
    checkpointer.loadFromCheckpoint(octDA,dtree,functionToRetain,vecs,domain, nullptr, false);
    IO::writeBoundaryElements(octDA, dtree.getTreePartFiltered(), "boundary_checkpoint", "boundary", domain);
    if(octDA->isActive()) {
      PETSc::petscVectopvtu(octDA, dtree.getTreePartFiltered(), vecs[0].v, "checkpoint", varname, domain, false, false,
                            ndof);
    }
  }
  else {
    TALYFEMLIB::PrintStatus("Storing checkpoint");
    octDA = createSubDA(dtree, functionToRetain, 5, eleOrder);
    subDA.finalize(octDA, dtree.getTreePartFiltered(), domain);
    for (int i = 0; i < 2; i++) {
      SDARefine sdaRefine(octDA, dtree.getTreePartFiltered(), domain, 20);
      DA *newDA = sdaRefine.getRefineSubDA(dtree);
      std::swap(newDA, octDA);
      delete newDA;
      subDA.finalize(octDA, dtree.getTreePartFiltered(), domain);
    }


    Vec tempVec;
    octDA->petscCreateVector(tempVec, false, false, ndof);
    std::function<void(const double *, double *)> initial_condition = [](const double *x, double *var) {
      var[0] = x[0]+x[1];
      var[1] = 20;
    };
    if(octDA->isActive()) {
      octDA->petscSetVectorByFunction(tempVec, initial_condition, false, false, ndof);
    }
    Checkpointer checkpointer(5, "Checkpoint");

    VecInfo Vec(tempVec, ndof, 0);
    std::vector<VecInfo> vecs(1, Vec);
    checkpointer.storeCheckpoint(octDA, &dtree, vecs, domain);
  if(octDA->isActive()) {
    PETSc::petscVectopvtu(octDA, dtree.getTreePartFiltered(), tempVec, "original", varname, domain, false, false, ndof);
  }

  }
  dendrite_finalize(octDA);


}