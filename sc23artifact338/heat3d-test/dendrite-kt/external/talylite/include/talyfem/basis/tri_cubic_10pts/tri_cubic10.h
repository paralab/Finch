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
#include <talyfem/basis/tri_cubic_10pts/tri_cubic10_impl.h>

namespace TALYFEMLIB {

/**
 * Cubic basis function for 2D triangle elements (10 points).
 */
template <int surface_id = 0, int rel_order = 0>
struct TriCubic10Basis : public BasisCommon< TriCubic10BasisImpl, TriItgPts<4 + rel_order, surface_id> > {
};

}  // namespace TALYFEMLIB
