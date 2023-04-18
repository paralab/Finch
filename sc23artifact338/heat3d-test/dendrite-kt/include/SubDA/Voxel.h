//
// Created by maksbh on 7/18/20.
//

#ifndef DENDRITEKT_VOXEL_H
#define DENDRITEKT_VOXEL_H

#include <DataTypes.h>
#include <Geometry/Geometry.h>
namespace VOXEL{

struct Sphere{
  DENDRITE_REAL center[3];
  DENDRITE_REAL radius;
  const RetainSide retainSide;
  Sphere(const DENDRITE_REAL _center[3], const DENDRITE_REAL _radius, const RetainSide _retainSide = RetainSide::OUT)
  :retainSide(_retainSide){
    std::memcpy(center,_center,sizeof(DENDRITE_REAL)*3);
    radius = _radius;
  }

  bool ifInside(const DENDRITE_REAL * pos) const{
    double xDist = pos[0] - center[0];
    double yDist = pos[1] - center[1];
    double zDist = pos[2] - center[2];
    double dist = xDist * xDist + yDist * yDist + zDist * zDist;
    bool isInside = false;
    if (dist < radius * radius) {
      isInside = true;
    }
    isInside = (retainSide==RetainSide::OUT) ? isInside:not(isInside);
    return isInside;
  }
};

struct Circle{
  DENDRITE_REAL center[2];
  DENDRITE_REAL radius;
  const RetainSide retainSide;
  Circle(const DENDRITE_REAL _center[2],const  DENDRITE_REAL _radius,const RetainSide _retainSide = RetainSide::OUT)
  :retainSide(_retainSide){
    std::memcpy(center,_center,sizeof(DENDRITE_REAL)*2);
    radius = _radius;
  }

  bool ifInside(const DENDRITE_REAL * pos) const{
    double xDist = pos[0] - center[0];
    double yDist = pos[1] - center[1];
    double dist = xDist * xDist + yDist * yDist;
    bool isInside = false;
    if (dist < radius * radius) {
      isInside = true;
    }
    isInside = (retainSide==RetainSide::OUT) ? isInside:not(isInside);
    return isInside;
  }
};

struct Box{
  DENDRITE_REAL min[DIM];
  DENDRITE_REAL max[DIM];
  const RetainSide retainSide;
  Box(const DENDRITE_REAL _min[DIM],const DENDRITE_REAL _max[DIM],const RetainSide _retainSide = RetainSide::OUT)
  :retainSide(_retainSide){
    std::memcpy(min,_min,sizeof(DENDRITE_REAL)*DIM);
    std::memcpy(max,_max,sizeof(DENDRITE_REAL)*DIM);
  }
  bool ifInside(const DENDRITE_REAL * pos) const {
    bool isInside = true;
    for (int dim = 0; dim < DIM; dim++) {
      if ((pos[dim] < min[dim]) or (pos[dim] > max[dim])) {
        isInside = false;
      }
    }
    isInside = (retainSide==RetainSide::OUT) ? isInside:not(isInside);
    return isInside;
  }

};


class Voxel{
 protected:

  std::vector<Sphere> m_sphere;
  std::vector<Circle> m_circle;
  std::vector<Box> m_box;
  std::vector<const GEOMETRY::Geometry*> m_geometry;


 public:

  void addObject(const VOXEL::Sphere & sphere);
  void addObject(const VOXEL::Circle & circle);
  void addObject(const VOXEL::Box & box);
  void addObject(const GEOMETRY::Geometry * geometry);



  inline const std::vector<Sphere> & getVoxelSpheres() const{
    return m_sphere;
  };

  inline const std::vector<Circle> & getVoxelCircles() const{
    return m_circle;
  };

  inline const std::vector<Box> & getVoxelBox() const{
    return m_box;
  };

  inline const std::vector<const GEOMETRY::Geometry*> & getGeometry() const{
    return m_geometry;
  };


};
}
#endif //DENDRITEKT_VOXEL_H
