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
 * Implementation for quadratic basis function for triangle elements.
 */
class TriQuadraticBasisImpl {
 public:
  static constexpr int nsd = 2;  ///< number of spatial dimensions
  static constexpr int nbf = 6;  ///< number of shape functions
  static constexpr int nbf_per_node = 1;  ///< shape functions per node
  static constexpr double jacobian_scale = 0.5;   ///< scale jacobian by 1/2

  /**
   * Calculate shape functions.
   * @param localPt integration point
   * @param[out] n_out shape functions evaluated at localPt (output)
   */
  static void calc_N(const ZEROPTV& localPt, double (&n_out)[nbf]) {
    n_out[0] = 2 * (1 - localPt.x() - localPt.y()) * (0.5 - localPt.x() - localPt.y());
    n_out[1] = 2 * localPt.x() * (localPt.x() - 0.5);
    n_out[2] = 2 * localPt.y() * (localPt.y() - 0.5);
    n_out[3] = 4 * (1 - localPt.x() - localPt.y()) * localPt.x();
    n_out[4] = 4 * localPt.x() * localPt.y();
    n_out[5] = 4 * (1 - localPt.x() - localPt.y()) * localPt.y();
  }

  /**
   * Calculate derivative of N in isoparametric space.
   * @param localPt integration point
   * @param[out] dnde_out derivative of N in isoparametric space (output)
   */
  static void calc_dNde(const ZEROPTV& localPt, double (&dnde_out)[nbf][nsd]) {
    double x = localPt.x();
    double y = localPt.y();
    dnde_out[0][0] = 4 * x + 4 * y - 3;
    dnde_out[0][1] = 4 * x + 4 * y - 3;
    dnde_out[1][0] = 4 * x - 1;
    dnde_out[1][1] = 0;
    dnde_out[2][0] = 0;
    dnde_out[2][1] = 4 * y - 1;
    dnde_out[3][0] = 4 - 4 * y - 8 * x;
    dnde_out[3][1] = -4 * x;
    dnde_out[4][0] = 4 * y;
    dnde_out[4][1] = 4 * x;
    dnde_out[5][0] = -4 * y;
    dnde_out[5][1] = 4 - 8 * y - 4 * x;
  }

  /**
     * Calculate SECOND derivative of N in isoparametric space.
     * @param localPt integration point
     * @param[out] dnde_out derivative of N in isoparametric space (output)
     */
  /*static void calc_dN2de(const ZEROPTV& localPt, double (&dn2de_out)[nbf][nsd]) {
        dn2de_out[0][0] = 4;
        dn2de_out[0][1] = 0;
        dn2de_out[1][0] = 0; 
        dn2de_out[1][1] = 4;
        dn2de_out[2][0] = 4; 
        dn2de_out[2][1] = 4;
        dn2de_out[3][0] = 0; 
        dn2de_out[3][1] = 0;
        dn2de_out[4][0] = 0; 
        dn2de_out[4][1] =-8;
        dn2de_out[5][0] =-8; 
        dn2de_out[5][1] = 0;
  }*/
};

}  // namespace TALYFEMLIB
