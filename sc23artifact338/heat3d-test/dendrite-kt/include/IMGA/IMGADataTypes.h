//
// Created by maksbh on 2/26/21.
//

#ifndef DENDRITEKT_IMGADATATYPES_H
#define DENDRITEKT_IMGADATATYPES_H

#include <DataTypes.h>

enum ElementMarker : DENDRITE_UINT {
  IN_ELEMENT = 0,
  OUT_ELEMENT = 1,
  INTERCEPTED_ELEMENT = 2,
  IN_GP = 3,
  OUT_GP = 4,
  INTERCEPTED_GP= 5,
};
enum MarkerType: bool{
  ELEMENT_NODES = false,
  GAUSS_POINT = true
};

enum IBM_METHOD : DENDRITE_UINT {
  SBM = 0,
  NITSCHE = 1
};

#endif //DENDRITEKT_IMGADATATYPES_H
