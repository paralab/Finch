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

namespace TALYFEMLIB {

/**
 * Integration points for 1D line, 2D box, and 3D hexahedron.
 * Integration points are defined for a normalized isoparametric element.
 * Actual values are provided via template specialization in the files
 * included below.
 * @param order The "N" value, or the number of 1D points.
 *              Order is usually chosen by the order of the basis function
 *              (2 for linear, 3 for quadratic, 4 for cubic...).
 *              The total number of gauss points is order^nsd.
 * @param nsd Dimension of points (can be 1D/2D/3D)
 * @param surface_id Surface ID, as defined by ELEM::GetSurfaceCheckArray().
 *                   0 means volume integration. Other numbers pick a surface.
 *                   (-1 = X-, +1 = X+, -2 = Y-, +2 = Y+, -3 = Z-, +3 = Z+)
 *                   Surfaces are only valid for matching dimensions
 *                   (i.e. using -3 in 2D will give a compiler error).
 */
template <int order, int nsd, int surface_id>
struct BoxItgPts {
};

/**
 * Integration points for 2D triangle.
 * Integration points are defined for a normalized isoparametric element.
 * Actual values are provided via template specialization.
 * @param order "N" value - 2 for linear, 3 for quadratic, etc.
 * @param surface_id Surface ID, as defined by ELEM::GetSurfaceCheckArray().
 *                   0 means volume integration. Other numbers pick a surface.
 *                   (-1 = X-, -2 = Y-, +1 = hypotenuse)
 */
template <int order, int surface_id>
struct TriItgPts {
};

/**
 * Integration points for 3D tetrahedron.
 * Integration points are defined for a normalized isoparametric element.
 * Actual values are provided via template specialization.
 * @param order "N" value - 2 for linear, 3 for quadratic, etc.
 * @param surface_id Surface ID, as defined by ELEM::GetSurfaceCheckArray().
 *                   0 means volume integration. Other numbers pick a surface.
 *                   (1 = hypotenuse, 2 = X-, 3 = Y-, 4 = Z-)
 */
template <int order, int surface_id>
struct TetItgPts {
};


#ifdef ENABLE_4D
/**
 * Integration points for 4D pentatope.
 * Integration points are defined for a normalized isoparametric element.
 * Actual values are provided via template specialization.
 * @param order "N" value - 2 for linear, 3 for quadratic, etc.
 * @param surface_id Surface ID, as defined by ELEM::GetSurfaceCheckArray().
 *                   0 means volume integration. Other numbers pick a surface.
 *                   (1 = hypotenuse, 2 = X-, 3 = Y-, 4 = Z-)
 */
template <int order, int surface_id>
struct PentItgPts {
};
#endif

}  // namespace TALYFEMLIB

#include <talyfem/basis/itg_pts/box_1d.h>
#include <talyfem/basis/itg_pts/box_2d.h>
#include <talyfem/basis/itg_pts/box_3d.h>
#include <talyfem/basis/itg_pts/box_4d.h>
#include <talyfem/basis/itg_pts/box_surfaces.h>

#include <talyfem/basis/itg_pts/tri_2d.h>
#include <talyfem/basis/itg_pts/tet.h>

#ifdef ENABLE_4D
#include <talyfem/basis/itg_pts/pent.h>
#endif