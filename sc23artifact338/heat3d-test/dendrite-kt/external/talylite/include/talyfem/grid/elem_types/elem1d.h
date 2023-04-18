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
#ifndef GRID_ELEM_TYPES_ELEM1D_H_
#define GRID_ELEM_TYPES_ELEM1D_H_

#include <talyfem/grid/elem.h>  // parent class


namespace TALYFEMLIB {

/**
 * One dimensional element.
 */
class ELEM1d : public ELEM {
 public:
  ElemType elmType() const override { return kElem1d; }
  int GetNodesPerSurface() const override;
  const int* GetSurfaceCheckArray() const override;
  int GetSurfaceCount() const override { return 2; }
  int nsd() const override { return 1; }
  ZEROPTV CalculateNormal(const GRID* grid, int surface_id) const override;
  double GetMeasure(const GRID* grid) const override;
  kBasisFunction basis_function() const override;
};

}  // namespace TALYFEMLIB

#endif  // GRID_ELEM_TYPES_ELEM1D_H_
