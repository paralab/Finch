////
//// Created by maksbh on 7/13/20.
////
//
//#ifndef DENDRITEKT_GEOMTYPE_H
//#define DENDRITEKT_GEOMTYPE_H
//
//#include <DataTypes.h>
//
//class GeomType{
// public:
//#if(DIM == 3)
// enum Type:DENDRITE_UINT {
//   OTHERS = 0,
//   SPHERE = 1
//  };
//
// struct Sphere{
//   const DENDRITE_REAL radius;
//   const DENDRITE_REAL center[3];
//
//   /**
//    * @brief constructor
//    * @param _radius radius of the sphere
//    * @param _center center of the sphere
//    */
//   Sphere(const DENDRITE_REAL _radius, const DENDRITE_REAL _center[3])
//   :radius(_radius),center{_center[0],_center[1],_center[2]}{
//
//   }
//   /**
//    * @param position the 3D position
//    * @param [out] normal returns the analytically calculated outward normal (must be allocated outside)
//    */
//   void getNormal(const DENDRITE_REAL * position, DENDRITE_REAL *normal) const{
//     /// TODO : Check it. Not sure if its correct
//     normal[0] = position[0] - center[0];
//     normal[1] = position[1] - center[1];
//     normal[2] = position[2] - center[2];
//     DENDRITE_REAL  mag = normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2];
//     normal[0] /= mag;
//     normal[1] /= mag;
//     normal[2] /= mag;
//   }
// };
//#endif
//};
//
//#endif //DENDRITEKT_GEOMTYPE_H
