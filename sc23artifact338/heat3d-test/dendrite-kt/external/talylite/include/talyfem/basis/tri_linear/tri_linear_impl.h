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
 * Implementation for linear basis function for triangle elements.
 */
class TriLinearBasisImpl {
 public:
  static constexpr int nsd = 2;  ///< number of spatial dimensions
  static constexpr int nbf = 3;  ///< number of shape functions
  static constexpr int nbf_per_node = 1;  ///< shape functions per node
  static constexpr double jacobian_scale = 0.5;  ///< scale jacobian by 1/2

  /**
   * Calculate shape functions.
   * @param localPt integration point
   * @param[out] n_out shape functions evaluated at localPt (output)
   */
  static void calc_N(const ZEROPTV& localPt, double (&n_out)[nbf]) {
    n_out[0] = 1 - localPt.x() - localPt.y();
    n_out[1] = localPt.x();
    n_out[2] = localPt.y();
  }

  /**
   * Calculate derivative of N in isoparametric space.
   * @param localPt integration point
   * @param[out] dnde_out derivative of N in isoparametric space (output)
   */
  static void calc_dNde(const ZEROPTV& localPt, double (&dnde_out)[nbf][nsd]) {
    dnde_out[0][0] = -1.0;
    dnde_out[0][1] = -1.0;

    dnde_out[1][0] = 1.0;
    dnde_out[1][1] = 0.0;

    dnde_out[2][0] = 0.0;
    dnde_out[2][1] = 1.0;
  }
};

}  // namespace TALYFEMLIB
