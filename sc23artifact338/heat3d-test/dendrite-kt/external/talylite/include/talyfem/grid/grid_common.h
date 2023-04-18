/*
  Copyright 2014-2016 Baskar Ganapathysubramanian

  This file is part of TALYFem.

  TALYFem is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as
  published by the Free Software Foundation, either version 2.1 of the
  License, or (at your option) any later version.

  TALYFem is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with TALYFem.  If not, see <http://www.gnu.org/licenses/>.
*/
// --- end license text --- //
#ifndef GRID_GRID_COMMON_H_
#define GRID_GRID_COMMON_H_

#include <talyfem/grid/elem_common.h>


namespace TALYFEMLIB {

/**
 * Method of handling parallel processing
 */
enum ParallelMethod {
  kNoDomainDecomp = 1,  ///< does not use domain decomposition
  kWithDomainDecomp = 2,  ///< uses domain decomposition
};


/**
 * Type of Grid
 */
enum GridType {
  kGrid3dBox = 0,  ///< 3D Box Grid
  kGrid2dBox = 1,  ///< 2D Box Grid
  kGrid2dTriangle = 2,  ///< 2D Triangle Grid
  kGrid3dTet = 3,  ///< 3D Tet Grid
  kGrid1d = 4,  ///< 1D Grid
};

/**
 * @param t grid type
 * @returns the associated element type for the given grid type
 */
inline ElemType grid_to_elem_type(GridType t) {
  switch (t) {
    case kGrid3dBox:
      return kElem3dHexahedral;
    case kGrid2dBox:
      return kElem2dBox;
    case kGrid2dTriangle:
      return kElem2dTriangle;
    case kGrid3dTet:
      return kElem3dTetrahedral;
    case kGrid1d:
      return kElem1d;
  }

  throw TALYException() << "Unknown GridType " << t << "!";
}

/**
 * @param t element type
 * @returns the grid type for elements of type t
 */
inline GridType elem_to_grid_type(ElemType t) {
  switch (t) {
    case kElem3dHexahedral:
      return kGrid3dBox;
    case kElem2dBox:
      return kGrid2dBox;
    case kElem2dTriangle:
      return kGrid2dTriangle;
    case kElem3dTetrahedral:
      return kGrid3dTet;
    case kElem1d:
      return kGrid1d;
  }

  throw TALYException() << "Unknown GridType " << t << "!";
}

}  // namespace TALYFEMLIB

#endif  // GRID_GRID_COMMON_H_
