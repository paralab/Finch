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
#pragma once

#include <talyfem/basis/basis_common.h>
#include <talyfem/basis/itg_pts/itg_pts.h>
#include <talyfem/basis/tri_quadratic/tri_quadratic_impl.h>

namespace TALYFEMLIB {

/**
 * Quadratic basis function for 2D triangle elements.
 */
template <int surface_id = 0, int rel_order = 0>
struct TriQuadraticBasis : public BasisCommon< TriQuadraticBasisImpl, TriItgPts<3 + rel_order, surface_id> > {
};

}  // namespace TALYFEMLIB
