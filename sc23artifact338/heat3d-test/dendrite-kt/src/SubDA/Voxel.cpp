//
// Created by maksbh on 7/18/20.
//

#include <SubDA/Voxel.h>
#include <DataTypes.h>

namespace VOXEL {

void Voxel::addObject(const VOXEL::Sphere &sphere) {
#if(DIM == 3)
  m_sphere.push_back(sphere);
#endif
}

void Voxel::addObject(const VOXEL::Circle &circle) {
#if(DIM == 2)
  m_circle.push_back(circle);
#endif
}

void Voxel::addObject(const VOXEL::Box &box) {
  m_box.push_back(box);
}

void Voxel::addObject(const GEOMETRY::Geometry * geometry){
//#if(DIM == 3)
  m_geometry.push_back(geometry);
//#endif

}
}
