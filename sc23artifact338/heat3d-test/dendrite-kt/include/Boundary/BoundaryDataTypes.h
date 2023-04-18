//
// Created by maksbh on 7/23/20.
//

#ifndef DENDRITEKT_BOUNDARYDATATYPES_H
#define DENDRITEKT_BOUNDARYDATATYPES_H
#include <DataTypes.h>
namespace BoundaryTypes {
enum WALL : DENDRITE_UINT {
    X_MINUS = 0,
    X_PLUS = 1,
    Y_MINUS = 2,
    Y_PLUS = 3,

#if (DIM == 3)
    Z_MINUS = 4,
    Z_PLUS =  5,
#endif
  MAX_WALL_TYPE_BOUNDARY = 6,
};

enum VOXEL : DENDRITE_UINT {
  SPHERE = MAX_WALL_TYPE_BOUNDARY + 0,
  BOX =    MAX_WALL_TYPE_BOUNDARY + 1,
  CIRCLE = MAX_WALL_TYPE_BOUNDARY + 2,
  GEOMETRY = MAX_WALL_TYPE_BOUNDARY + 3,
  MAX_BOUNDARY_TYPE = MAX_WALL_TYPE_BOUNDARY + 4,
};
}
#define MAX_BOUNDARY_TYPES 12
static_assert(BoundaryTypes::VOXEL::MAX_BOUNDARY_TYPE <= MAX_BOUNDARY_TYPES,"Error in max boundary types");


#endif //DENDRITEKT_BOUNDARYDATATYPES_H
