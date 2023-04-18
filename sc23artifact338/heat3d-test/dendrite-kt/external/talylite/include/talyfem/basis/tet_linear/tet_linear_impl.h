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

#include <talyfem/grid/zeroptv.h>

namespace TALYFEMLIB {

/**
 * Implementation for linear basis function for tetrahedral elements.
 */
class TetLinearBasisImpl {
 public:
  static constexpr int nsd = 3;  ///< number of spatial dimensions
  static constexpr int nbf = 4;  ///< number of shape functions
  static constexpr int nbf_per_node = 1;  ///< shape funcs per node
  static constexpr double jacobian_scale = 1.0 / 6.0;  ///< scale jacobian by

  /**
   * Calculate shape functions.
   * @param localPt integration point
   * @param[out] n_out shape functions evaluated at localPt (output)
   */
  static void calc_N(const ZEROPTV& localPt, double (&n_out)[nbf]) {
    n_out[0] = 1 - localPt.x() - localPt.y() - localPt.z();
    n_out[1] = localPt.x();
    n_out[2] = localPt.y();
    n_out[3] = localPt.z();
  }

  /**
   * Calculate derivative of N in isoparametric space.
   * @param localPt integration point
   * @param[out] dnde_out derivative of N in isoparametric space (output)
   */
  static void calc_dNde(const ZEROPTV& localPt, double (&dnde_out)[nbf][nsd]) {
    dnde_out[0][0] = -1.0;
    dnde_out[0][1] = -1.0;
    dnde_out[0][2] = -1.0;

    dnde_out[1][0] = 1.0;
    dnde_out[1][1] = 0.0;
    dnde_out[1][2] = 0.0;

    dnde_out[2][0] = 0.0;
    dnde_out[2][1] = 1.0;
    dnde_out[2][2] = 0.0;

    dnde_out[3][0] = 0.0;
    dnde_out[3][1] = 0.0;
    dnde_out[3][2] = 1.0;
  }
};

}  // namespace TALYFEMLIB
