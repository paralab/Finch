//
// Created by maksbh on 7/18/20.
//

#ifndef DENDRITEKT_SDAREFINE_H
#define DENDRITEKT_SDAREFINE_H

#include <Traversal/Refinement.h>
#include <SubDA/Voxel.h>

class SDARefine : public Refinement {

  const DENDRITE_UINT blockRefinement_;
  const DENDRITE_UINT objectRefinement_;
  const VOXEL::Voxel & voxel_;
 public:
  SDARefine(DA *da,  const std::vector<TREENODE> & treePart, const DomainExtents &domainInfo, const VOXEL::Voxel &voxel,
            const DENDRITE_UINT blockRefinement, const double objectRefinement);

  DENDRITE_UINT calcSphereRefineLevel(const std::vector<TALYFEMLIB::ZEROPTV> &coords);
  DENDRITE_UINT calcCircleRefineLevel(const std::vector<TALYFEMLIB::ZEROPTV> &coords);
  DENDRITE_UINT calcBoxRefineLevel(const std::vector<TALYFEMLIB::ZEROPTV> &coords);

  ot::OCT_FLAGS::Refine getRefineFlags(TALYFEMLIB::FEMElm &fe, const std::vector<TALYFEMLIB::ZEROPTV> &coords) override;


};

SDARefine::SDARefine(DA *da, const std::vector<TREENODE> & treePart,const DomainExtents &domainInfo,  const VOXEL::Voxel &voxel,
                     const DENDRITE_UINT blockRefinement, const double objectRefinement)
    : Refinement(da, treePart,domainInfo), blockRefinement_(blockRefinement),
      objectRefinement_(objectRefinement), voxel_(voxel) {

  this->traverse();
}

DENDRITE_UINT SDARefine::calcBoxRefineLevel(const std::vector<TALYFEMLIB::ZEROPTV> &coords) {
  DENDRITE_UINT reqLevelBox = this->m_level;
  const std::vector<VOXEL::Box> & box = voxel_.getVoxelBox();
  double boxPadWidth = 0.5;
  bool foundForRefinement = false;
  for (int i = 0; i < box.size(); i++) {
    for (int numCoords = 0; numCoords < 8; numCoords++) {
      if ((coords[numCoords][0] > (box[i].min[0] - boxPadWidth)) and (coords[numCoords][0] < (box[i].max[0] + boxPadWidth)) and
          ((coords[numCoords][1] > (box[i].min[1] - boxPadWidth)) and (coords[numCoords][1] < (box[i].max[1] + boxPadWidth))))
      {
        return objectRefinement_;
      }


    }
  }
  return reqLevelBox;
}

DENDRITE_UINT SDARefine::calcSphereRefineLevel(const std::vector<TALYFEMLIB::ZEROPTV> &coords) {
  DENDRITE_UINT reqLevelSphere = this->m_level;
  bool foundForRefinement = false;
  const std::vector<VOXEL::Sphere> & spheres = voxel_.getVoxelSpheres();
  for (int i = 0; i < spheres.size(); i++) {
    if (not(foundForRefinement)) {
      const double spherePadWidth = 1.0 * spheres[i].radius; /// This needs to be properly optimized
      for (int numCoords = 0; numCoords < 8; numCoords++) {
        double distX = spheres[i].center[0] - coords[numCoords].x();
        double distY = spheres[i].center[1] - coords[numCoords].y();
        double distZ = spheres[i].center[2] - coords[numCoords].z();
        double dist = distX * distX + distY * distY + distZ * distZ;
        if (dist < spherePadWidth * spherePadWidth) {
          foundForRefinement = true;
          reqLevelSphere = objectRefinement_;
          break;
        }
      }
    }
    if(foundForRefinement){
      break;
    }
  }
  return reqLevelSphere;
}

DENDRITE_UINT SDARefine::calcCircleRefineLevel(const std::vector<TALYFEMLIB::ZEROPTV> &coords) {
  DENDRITE_UINT reqLevelCircle = this->m_level;
  bool foundForRefinement = false;
  const std::vector<VOXEL::Circle> & circles = voxel_.getVoxelCircles();
  for (int i = 0; i < circles.size(); i++) {
    if (not(foundForRefinement)) {
      const double circlePadWidth = 1.1 * circles[i].radius; /// This needs to be properly optimized
      for (int numCoords = 0; numCoords < 4; numCoords++) {
        double distX = circles[i].center[0] - coords[numCoords].x();
        double distY = circles[i].center[1] - coords[numCoords].y();
        double dist = distX * distX + distY * distY;
        if (dist < circlePadWidth*circlePadWidth) {
          foundForRefinement = true;
          reqLevelCircle = objectRefinement_;
          break;
        }
      }
    }
    if(foundForRefinement){
      break;
    }
  }
  return reqLevelCircle;
}

ot::OCT_FLAGS::Refine
SDARefine::getRefineFlags(TALYFEMLIB::FEMElm &fe, const std::vector<TALYFEMLIB::ZEROPTV> &coords) {
//  return ot::OCT_FLAGS::Refine::OCT_REFINE;
  enum RefineLevel{
    BASE = 0,
    BLOCK = 1,
    CIRCLE = 2,
    SPHERE = 3,
    BOX = 4,
    MAX = 5
  };
  std::vector<DENDRITE_UINT> refineLevel(RefineLevel::MAX,this->m_level);

  DENDRITE_UINT currLevel = this->m_level;

  if (coords[7].z() >= 140 and coords[0].z() <= 170) {
    refineLevel[RefineLevel::BLOCK] = blockRefinement_;
  }

  refineLevel[RefineLevel::SPHERE] = calcSphereRefineLevel(coords);
  refineLevel[RefineLevel::CIRCLE] = calcCircleRefineLevel(coords);
  refineLevel[RefineLevel::BOX] = calcBoxRefineLevel(coords);
  DENDRITE_UINT reqLevel = *std::max_element(refineLevel.begin(),refineLevel.end());


  if (this->m_BoundaryOctant) {
    return ot::OCT_FLAGS::Refine::OCT_REFINE;
  }
  return ot::OCT_FLAGS::Refine::OCT_NO_CHANGE;
}


#endif //DENDRITEKT_SDAREFINE_H
