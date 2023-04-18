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

#include <talyfem/basis/basis.h>
#include <talyfem/basis/elemnodes.h>
#include <talyfem/grid/zeroptv.h>

namespace TALYFEMLIB {

/**
 * Basis function for "points" (1D surfaces).
 */
class PointBasis {
 public:
  static const int nsd = 0;  ///< number of spatial dimensions
  static const int nbf = 1;  ///< number of basis functions
  static const int nbf_per_node = 1;  ///< bf per node
  static const int n_itg_pts = 1;  ///< number of itg points

  /**
   * @param i ignored
   * @returns ZEROPTV(0, 0, 0)
   */
  static constexpr ZEROPTV itg_pts(int i) { return ZEROPTV(0, 0, 0); }

  /**
   * @param i gauss point index (ignored)
   * @returns 1
   */
  static constexpr double weights(int i) { return 1.0; }

  /**
   * @param itg_pt itg pt index
   * @param elem node positions
   * @param flags ignored
   * @param[out] vals output
   * @param[inout] rotation_matrix ignored
   * @param[inout] elem_cache ignored
   */
  static void calc_values(int itg_pt, const ElemNodes& elem,
                          unsigned int flags, BasisValues<nbf, nsd>* vals,
                          Matrix3& rotation_matrix, ElemNodes& elem_cache) {
    vals->itg_pt = itg_pts(0);
    vals->N[0] = 1.0;
    vals->jacobian = 1.0;
    vals->jacc_x_weight = 1.0;
    vals->position = elem.node_pos(0);
  }

  /**
   * @param itg_pt itg pt index
   * @param elem node positions
   * @param flags ignored
   * @param[out] vals output
   */
  static void calc_values(const ZEROPTV& itg_pt, const ElemNodes& elem,
                          unsigned int flags, BasisValues<nbf, nsd>* vals) {
    Matrix3 dummy;
    ElemNodes dummy2;
    return calc_values(0, elem, flags, vals, dummy, dummy2);
  }
};

}  // namespace TALYFEMLIB
