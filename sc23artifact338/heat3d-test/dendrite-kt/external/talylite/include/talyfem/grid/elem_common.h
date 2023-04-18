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
#ifndef GRID_ELEM_COMMON_H_
#define GRID_ELEM_COMMON_H_

#include <talyfem/utils/utils.h>  // for PrintError
#include <talyfem/common/exceptions.h>

//#ifdef ENABLE_4D
//#include <talyfem/grid/elem_types/elem4dtesseract.h>
//#endif


namespace TALYFEMLIB {

/**
 * Type of finite element
 */
enum ElemType {
  kElem3dHexahedral = 0,  // 3D Box Grid
  kElem2dBox = 1,  // 2D Box Grid
  kElem2dTriangle = 2,  // 2D Triangle Grid
  kElem3dTetrahedral = 3,  // 3D Tet Grid
  kElem1d = 4,  // 1D Grid
#ifdef ENABLE_4D
  kElem4dTesseract = 6, // 4D Tesseract Grid
#endif
};

/**
 * Get the number of nodes per element.
 * Depends on the order of N in the grid.
 *
 * @param elem_type Element type.
 * @param order Order of N, e.g. grid->getOrderOfN().
 * @returns number of nodes in the given element.
 */
inline unsigned int get_nodes_in_element(ElemType elem_type, int order = 1) {
  switch (elem_type) {
    case kElem3dHexahedral:
      return static_cast<unsigned int>((order + 1) * (order + 1) * (order + 1));
    case kElem2dBox:
      return static_cast<unsigned int>((order + 1) * (order + 1));
    case kElem2dTriangle:
      if (order != 1) {
        throw NotImplementedException() <<
            "Only order 1 is supported for 2D triangle!";
      }
      return 3;  // only order 1 is supported
    case kElem3dTetrahedral:
      if (order != 1) {
        throw NotImplementedException() <<
            "Only order 1 is supported for 3D tetrahedral!";
      }
      return 4;  // only order 1 is supported
    case kElem1d:
      return static_cast<unsigned int>(1 + order);
#ifdef ENABLE_4D
    case kElem4dTesseract:
      return static_cast<unsigned int>((order + 1) * (order + 1) * (order + 1) * (order + 1));
#endif
  }

  return 0;
}

/**
 * Returns the number of spatial dimensions needed for the given element type.
 *
 * @param t ElemType to get information about
 * @return number of spatial dimensions for the specified element type
 */
inline unsigned int get_elem_type_nsd(ElemType t) {
  switch (t) {
    case kElem3dHexahedral: return 3;
    case kElem2dBox: return 2;
    case kElem2dTriangle: return 2;
    case kElem3dTetrahedral: return 3;
    case kElem1d: return 1;
#ifdef ENABLE_4D
    case kElem4dTesseract: return 4;
#endif
  }

  return 0;
}

}  // namespace TALYFEMLIB

#endif  // GRID_ELEM_COMMON_H_
