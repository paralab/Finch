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
#ifndef GRID_ELEM_TYPES_H_
#define GRID_ELEM_TYPES_H_

#include <talyfem/grid/elem_common.h>  // for ElemType
#include <talyfem/grid/elem_types/elem1d.h>  // for ELEM1d class
#include <talyfem/grid/elem_types/elem2dbox.h>  // for ELEM2dBox class
#include <talyfem/grid/elem_types/elem2dtriangle.h>  // for ELEM2dTriangle class
#include <talyfem/grid/elem_types/elem3dhexahedral.h>  // for ELEM3dHexahedral class
#include <talyfem/grid/elem_types/elem3dtetrahedral.h>  // for ELEM3dTetrahedral class

#ifdef ENABLE_4D
//#include <talyfem/grid/elem_types/elem4dpentatope.h>  // for ELEM4DPentatope class
  #include <talyfem/grid/elem_types/elem4dtesseract.h> // for ELEM4dTesseract class
#endif

namespace TALYFEMLIB {

class ELEM;

/**
 * @param type type of element to create
 * @returns a new element of the given type
 */
ELEM* make_elem_of_type(ElemType type);

}  // namespace TALYFEMLIB

#endif  // GRID_ELEM_TYPES_H_
