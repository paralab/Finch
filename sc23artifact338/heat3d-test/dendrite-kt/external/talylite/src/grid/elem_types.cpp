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
#include <talyfem/grid/elem-types.h>
#ifdef ENABLE_4D
#include <talyfem/grid/elem_types/elem4dtesseract.h>
#endif
namespace TALYFEMLIB {

ELEM* make_elem_of_type(ElemType type) {
  switch (type) {
    case kElem3dHexahedral:
      return new ELEM3dHexahedral();
    case kElem2dBox:
      return new ELEM2dBox();
    case kElem2dTriangle:
      return new ELEM2dTriangle();
    case kElem3dTetrahedral:
      return new ELEM3dTetrahedral();
    case kElem1d:
      return new ELEM1d();
#ifdef ENABLE_4D
    case kElem4dTesseract:
      return new ELEM4dTesseract();
#endif
  }

  PrintError("Unknown element type ", type);
  return NULL;
}

}  // namespace TALYFEMLIB
