

#include <PointCloud/PointData.h>
#include <DataTypes.h>
#include <DendriteUtils.h>
#include <PETSc/IO/petscVTU.h>
using namespace PETSc;
int main(int argc, char *argv[]) {
  dendrite_init(argc,argv);
  int eleOrder = 1;
  int level = 3;
  const DENDRITE_UINT  nx = 64;
  const DENDRITE_UINT  ny = 72;

  DomainInfo domainInfo;
  domainInfo.min.fill(0.0);
  domainInfo.max.fill(1.0);
  domainInfo.max[0] = 1.0;
  domainInfo.max[1] = 0.52;
  std::array<DENDRITE_REAL,DIM> dX{domainInfo.max[0]/(nx - 1.0), domainInfo.max[1]/(ny - 1.0)};

  DENDRITE_REAL position[DIM];
  PointCloud::PointData<PointCloud::InputType::CARTESIAN> pointData(dX,domainInfo,1);
  for(DENDRITE_UINT i = 0; i < nx; i++){
    for(DENDRITE_UINT j = 0; j < ny; j++){
      position[0] = domainInfo.min[0] + i*dX[0];
      position[1] = domainInfo.min[1] + j*dX[1];
      double val = position[0] + position[1];
      DENDRITE_UINT id[2]{i,j};
      pointData.addPoint(id,&val);
    }
  }


  std::vector<TREENODE> treeNode;
//  DA * octDA = createRegularDA(treeNode,level,eleOrder);
//  Vec vec;
//  setVectorFromPoints(octDA,pointData,domainInfo,vec,false,false);

  const char *varname[]{"u"};
//  petscVectopvtu(octDA,treeNode,vec,"fileOld",varname,domainInfo,false, false,1);

    DA * octDAnew = createRegularDA(treeNode,level,eleOrder);
    Vec vecnew;
    setVectorFromPoints(octDAnew,pointData,domainInfo,vecnew,false,false);
    pointData.finalize();
    petscVectopvtu(octDAnew,treeNode,vecnew,"fileNew",varname,domainInfo,false, false,1);
    dendrite_finalize(octDAnew);
}